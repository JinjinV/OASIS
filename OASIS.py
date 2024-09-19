 # -*- coding: utf-8 -*-

from osgeo import gdal, ogr, osr
import numpy as np
from collections import Counter
from sympy import *
from matplotlib import pyplot as plt
from os.path import exists
from os.path import basename
from os.path import splitext
import os, shutil, glob, subprocess, copy, sys, operator
import dbf
from adjustText import adjust_text
import time, pickle
from tqdm import tqdm

np.warnings.filterwarnings('ignore')

def raster2array(File):
    metadata = {}
    dataset = gdal.Open(File)
    metadata['array_rows'] = dataset.RasterYSize
    metadata['array_cols'] = dataset.RasterXSize
    metadata['geotransform'] = dataset.GetGeoTransform()
    metadata['projection'] = dataset.GetProjection()

    raster = dataset.GetRasterBand(1)
    metadata['noDataValue'] = raster.GetNoDataValue()
    metadata['scaleFactor'] = raster.GetScale()

    mapinfo = dataset.GetGeoTransform()
    metadata['pixelWidth'] = mapinfo[1]
    metadata['pixelHeight'] = mapinfo[5]

    metadata['bandstats'] = {}
    stats = raster.GetStatistics(True, True)
    metadata['bandstats']['min'] = round(stats[0], 2)
    metadata['bandstats']['max'] = round(stats[1], 2)
    metadata['bandstats']['mean'] = round(stats[2], 2)
    metadata['bandstats']['stdev'] = round(stats[3], 2)

    array = dataset.GetRasterBand(1).ReadAsArray(0, 0, metadata['array_cols'], metadata['array_rows']).astype(
        np.float)
    array[array == metadata['noDataValue']] = np.nan

    if not metadata['scaleFactor'] == None:
        array = array / metadata['scaleFactor']

    return array, metadata

def Directionlist(x, y, list, Raster):
    vallist = list
    if len(vallist) < 2:
        loop = 6
    else:
        loop = 2

    for k in range(0, loop):
        val = vallist[-1]
        if val == np.nan:
            break
        if val == 1:
            x += 1
            vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
        if val == 2:
            x += 1
            y -= 1
            vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
        if val == 3:
            y -= 1
            vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
        if val == 4:
            x -= 1
            y -= 1
            vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
        if val == 5:
            x -= 1
            vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
        if val == 6:
            x -= 1
            y += 1
            vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
        if val == 7:
            y += 1
            vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
        if val == 8:
            x += 1
            y += 1
            vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
    return vallist

def FindMajorDirection(D8raster, point):
    shp1 = ogr.Open(point)
    layer1 = shp1.GetLayer()
    for feat in layer1:
        point = str(feat.GetGeometryRef())
    X = float(point.split()[1][1:])
    Y = float(point.split()[2][:-1])

    src_ds = gdal.Open(D8raster)
    gt = src_ds.GetGeoTransform()
    Raster = src_ds.GetRasterBand(1)

    x = int((float(X) - gt[0]) / gt[1])
    y = int((float(Y) - gt[3]) / gt[5])

    vallist = []
    vallist.append(float(Raster.ReadAsArray(x, y, 1, 1)))
    vallist = Directionlist(x, y, vallist, Raster)

    valmocom = Counter(vallist).most_common()
    try:
        if valmocom[0][1] == valmocom[1][1]:
            while valmocom[0][1] == valmocom[1][1]:
                vallist = Directionlist(x, y, vallist)
                valmocom = Counter(vallist).most_common()
        RiverDirection = valmocom[0][0]

    except:
        RiverDirection = valmocom[0][0]

    return X, Y, RiverDirection

def createDS(ds_name, ds_format, geom_type, srs, overwrite=False):
    drv = ogr.GetDriverByName(ds_format)
    if os.path.exists(ds_name) and overwrite is True:
        deleteDS(ds_name)
    ds = drv.CreateDataSource(ds_name)
    lyr_name = os.path.splitext(os.path.basename(ds_name))[0]
    lyr = ds.CreateLayer(lyr_name, srs, geom_type)
    return ds, lyr

def polygonize(rasterTemp, outShp):

    sourceRaster = gdal.Open(rasterTemp)
    band = sourceRaster.GetRasterBand(1)
    driver = ogr.GetDriverByName("ESRI Shapefile")

    if os.path.exists(outShp):
        driver.DeleteDataSource(outShp)

    outDatasource = driver.CreateDataSource(outShp)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(sourceRaster.GetProjectionRef())
    outLayer = outDatasource.CreateLayer(outShp, srs)

    newField = ogr.FieldDefn('Class', ogr.OFTInteger)
    outLayer.CreateField(newField)

    gdal.Polygonize(band, None, outLayer, 0, [], callback=None)

    outDatasource.Destroy()

    sourceRaster = None
    band = None

    OShpFile = ogr.Open(outShp, update=1)

    lyr = OShpFile.GetLayerByIndex(0)

    lyr.ResetReading()

    for i in lyr:
        lyr.SetFeature(i)

        if i.GetField('Class') < 0:
            lyr.DeleteFeature(i.GetFID())
    OShpFile.Destroy()

def FindDNStreamline(treedat):
    EndStreamList = []
    with open(treedat) as Fopen:
        Flines = Fopen.read().splitlines()
        for k in range(0, len(Flines)):
            if int(Flines[k].split("\t")[5]) == -1 and int(Flines[k].split("\t")[6]) == -1:
                EndStreamList.append(int(Flines[k].split("\t")[1]))

    DNStreams = []
    for LN in EndStreamList:
        DNStreams.append(GetDNStreamNum(LN, treedat))
    return DNStreams

def extractValue(shp_in_file):

    '''
    Opens the input file, copies it into the output file, checking
    the filter.
    '''

    in_ds = ogr.Open(shp_in_file)
    if in_ds is None:
        print("Open failed.\n")
        sys.exit(1)
    in_lyr = in_ds.GetLayerByName(splitext(basename(shp_in_file))[0])
    if in_lyr is None:
        print("Error opening layer")
        sys.exit(1)
    driver_name = "ESRI Shapefile"
    drv = ogr.GetDriverByName(driver_name)
    if drv is None:
        print("%s driver not available.\n" % driver_name)
        sys.exit(1)
    lyr_def = in_lyr.GetLayerDefn()
    in_lyr.ResetReading()
    StreamDic = {}
    for feat in in_lyr:
        Linkno = int(feat.GetFieldAsString(feat.GetFieldIndex("LINKNO")))
        Length = feat.GetGeometryRef().Length()
        Drop = float(feat.GetFieldAsString(feat.GetFieldIndex("Drop")))
        StreamDic[Linkno] = [Length, Drop]
    in_ds = None
    return StreamDic

def Shed(tifpath, AD8, D8, pointI, bypoint=None):
    # MultiCoreNumber
    CREATE_NO_WINDOW = 0x08000000
    MCN = os.cpu_count()

    # Microsoft MPI - mpiexec.exe 경로
    MPIpath = r"C:\Program Files\Microsoft MPI\Bin"
    MultiRun = r"mpiexec.exe" + " -n " + str(MCN)

    # TauDEMpath
    TDpath = os.path.join(os.getcwd(), "TauDEM\\")

    # InputDEMpath
    OriDEM = tifpath
    Fname = OriDEM.split("\\")[-1].split(".")[0]

    # TempPath
    Temppath = os.path.join(os.getcwd(), "temp")

    # Mlp open
    os.system("cd " + MPIpath)

    # Threshhold
    temparray, tempmeta = raster2array(AD8)

    threshold = np.nanmax(temparray) * 0.005

    if bypoint != None:

        Fname += "_" + pointI

        AreaD8Tiff = os.path.join(Temppath, Fname + "_AreaD8.tif")

        subprocess.call(
            MultiRun + " " + TDpath + "AreaD8.exe -p " + D8 + " -o " + bypoint + " -nc -ad8 " + AreaD8Tiff, creationflags=CREATE_NO_WINDOW)

        plenTiff = os.path.join(Temppath, Fname + "_plen.tif")
        tlenTiff = os.path.join(Temppath, Fname + "_tlen.tif")
        gordTiff = os.path.join(Temppath, Fname + "_gord.tif")

        subprocess.call(
            MultiRun + " " + TDpath + "Gridnet.exe -p " + D8 + " -o " + bypoint + " -plen " + plenTiff + " -tlen "
            + tlenTiff + " -gord " + gordTiff, creationflags=CREATE_NO_WINDOW)

        srcTiff = os.path.join(Temppath, Fname + "_src.tif")

        subprocess.call(MultiRun + " " + TDpath + "Threshold.exe -ssa " + AreaD8Tiff + " -thresh " + str(
            threshold) + " -src " + srcTiff, creationflags=CREATE_NO_WINDOW)

        ordTiff = os.path.join(Temppath, Fname + "_ord.tif")
        treeTiff = os.path.join(Temppath, Fname + "_tree.dat")
        coordiff = os.path.join(Temppath, Fname + "_coord.dat")
        netTiff = os.path.join(Temppath, Fname + "_net.shp")
        wTiff = os.path.join(Temppath, Fname + "_w.tif")
        subprocess.call(
            MultiRun + " " + TDpath + "StreamNet.exe -fel " + tifpath + " -p " + D8 + " -ad8 " + AreaD8Tiff + " -src "
            + srcTiff + " -ord " + ordTiff + " -tree " + treeTiff + " -coord " + coordiff + " -net " + netTiff
            + " -w " + wTiff, creationflags=CREATE_NO_WINDOW)
        # print(MultiRun + " " + TDpath + "StreamNet.exe -fel " + tifpath + " -p " + D8 + " -ad8 " + AreaD8Tiff + " -src "
        #     + srcTiff + " -o " + bypoint + " -ord " + ordTiff + " -tree " + treeTiff + " -coord " + coordiff + " -net " + netTiff
        #     + " -w " + wTiff)
        polygonize(wTiff, os.path.join(Temppath, Fname + "_shed.shp"))
        # removelist = [ordTiff, coordiff, srcTiff, gordTiff, tlenTiff, plenTiff]
        #
        # for jj in removelist:
        #     os.remove(jj)
        if not os.path.exists(os.path.join(Temppath, Fname + "_net.prj")):
            shutil.copy(os.path.join(Temppath, Fname + "_shed.prj"), os.path.join(Temppath, Fname + "_net.prj"))

        return Fname

    else:
        print("plz input point shape file")

def GetUpperStreamNum(TreePath):
    filepath = TreePath
    UpperStreamNumlist = []
    OldTemplist = []
    with open(filepath) as Datfile:
        Data = Datfile.readlines()
        Datadic = {}
        for nline in range(0, len(Data)):
            NData = Data[nline].split("\t")
            Datadic[int(NData[1])] = [int(NData[5]), int(NData[6])]
            if int(NData[4]) == -1:
                TargetNum = int(NData[1])
    UpperStreamNumlist.append(Datadic[TargetNum][0])
    UpperStreamNumlist.append(Datadic[TargetNum][1])
    Templist = copy.deepcopy(UpperStreamNumlist)
    while len(Templist) != 0:
        for k in range(0, len(Templist)):
            if Datadic[Templist[k]][0] != -1:
                UpperStreamNumlist.append(Datadic[Templist[k]][0])
            if Datadic[Templist[k]][1] != -1:
                UpperStreamNumlist.append(Datadic[Templist[k]][1])
        OldTemplist += Templist
        Templist = copy.deepcopy(UpperStreamNumlist)
        for kk in range(0, len(OldTemplist)):
            Templist.remove(OldTemplist[kk])
    UpperStreamNumlist.append(TargetNum)
    return UpperStreamNumlist

def GetDNStreamNum(TargetNum, treedat):

    filepath = treedat
    DNStreamNumlist = []
    OldTemplist = []
    with open(filepath) as Datfile:
        Data = Datfile.readlines()
        Datadic = {}
        for nline in range(0, len(Data)):
            NData = Data[nline].split("\t")
            Datadic[int(NData[1])] = [int(NData[4])]
    DNStreamNumlist.append(Datadic[TargetNum][0])
    Templist = copy.deepcopy(DNStreamNumlist)
    if len(Templist) == 1 and Templist[0] == -1:
        DNStreamNumlist.remove(-1)
    else:
        while len(Templist) != 0:
            for k in range(0, len(Templist)):
                if Datadic[Templist[k]][0] != -1:
                    DNStreamNumlist.append(Datadic[Templist[k]][0])
            OldTemplist += Templist
            Templist = copy.deepcopy(DNStreamNumlist)
            for kk in range(0, len(OldTemplist)):
                Templist.remove(OldTemplist[kk])
    DNStreamNumlist.append(TargetNum)
    return DNStreamNumlist

def filter_func(value, filter_values):
    '''
    The custom filter function. In this case, we check that the value is in the
    value list, stripping the white spaces. In the case of numeric values, a
    comparaison could be done
    '''
    if int(value) in filter_values:
        return True
    else:
        return False

def warp(inshp, inraster):

    shp = ogr.Open(inshp)
    layer = shp.GetLayer()
    spatial_ref = layer.GetSpatialRef()

    union = ogr.Geometry(ogr.wkbPolygon)
    outraster = os.path.join(os.getcwd(), "temp\\temp.tif")

    for feat in layer:
        geom = feat.GetGeometryRef()
        union = union.Union(geom)

    newuni = union.Buffer(0)
    write_shapefile(newuni, inshp[:-4] + "T.shp", spatial_ref)
    OutTile = gdal.Warp(outraster, inraster, format="GTiff", multithread=True, cutlineDSName=inshp[:-4] + "T.shp",
                        cropToCutline=True, )
    OutTile = None
    Raster = None
    os.remove(inshp[:-4] + "T.shp")
    Cutarray, CutMeta = raster2array(outraster)
    os.remove(outraster)

    return Cutarray, CutMeta

def GetValueFromTif(tif, Shp):

    array, meta = raster2array(tif)

    Shp_filename = Shp

    gt = meta['geotransform']

    ds = ogr.Open(Shp_filename)
    lyr = ds.GetLayer()
    Valuelist = []
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom.ExportToWkt()[:10] == "MULTIPOINT":
            Totxt = geom.ExportToWkt()[12:-1]
        else:
            Totxt = geom.ExportToWkt()[7:-1]
        Pointlist = Totxt.split(",")

    for Point in Pointlist:
        mx, my = float(Point.split()[0]), float(Point.split()[1])


        px = int((mx - gt[0]) / gt[1])
        py = int((my - gt[3]) / gt[5])

        Valuelist.append(float(array[py, px]))

    return Valuelist

def GetValueFromRaster(array, meta, Shp):

    Shp_filename = Shp

    gt = meta['geotransform']

    ds = ogr.Open(Shp_filename)
    lyr = ds.GetLayer()
    Valuelist = []
    for feat in lyr:
        geom = feat.GetGeometryRef()
        if geom.ExportToWkt()[:10] == "MULTIPOINT":
            Totxt = geom.ExportToWkt()[12:-1]
        else:
            Totxt = geom.ExportToWkt()[7:-1]
        Pointlist = Totxt.split(",")

    for Point in Pointlist:
        mx, my = float(Point.split()[0]), float(Point.split()[1])


        px = int((mx - gt[0]) / gt[1])
        py = int((my - gt[3]) / gt[5])

        Valuelist.append(float(array[py, px]))

    return Valuelist

def intersects_L(shp1, shp2, outshp):
    outofrange = 0
    shp1 = ogr.Open(shp1)
    layer1 = shp1.GetLayer()

    shp2 = ogr.Open(shp2)
    layer2 = shp2.GetLayer()
    geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
    for feature in layer2:
        geomcol.AddGeometry(feature.GetGeometryRef())
    convexhull = geomcol.ConvexHull()

    union1 = ogr.Geometry(layer1.GetGeomType())
    for feat in layer1:
        geom = feat.GetGeometryRef()
        union1 = union1.Union(geom)

    intersection = union1.Intersection(convexhull)
    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(outshp)

    if intersection.ExportToWkt().split()[0] == "MULTIPOINT":
        layer = ds.CreateLayer('', None, ogr.wkbMultiPoint)
    elif intersection.ExportToWkt().split()[0] == "POINT":
        layer = ds.CreateLayer('', None, ogr.wkbPoint)
    elif intersection.ExportToWkt().split()[0] == "MULTILINESTRING":
        layer = ds.CreateLayer('', None, ogr.wkbMultiLineString)
    elif intersection.ExportToWkt().split()[0] == "LINESTRING":
        layer = ds.CreateLayer('', None, ogr.wkbLineString)
    else:
         outofrange = 1
    if outofrange != 1:
        layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        defn = layer.GetLayerDefn()
        feat = ogr.Feature(defn)
        feat.SetField('id', 123)
        geom = ogr.CreateGeometryFromWkt(intersection.ExportToWkt())
        feat.SetGeometry(geom)
        Length = geom.Length()
        layer.CreateFeature(feat)

        feat = geom = None

    ds = layer = feat = geom = None
    return Length

def intersects(Searching_Points, shp1, shp2, Buffdistance):
    # Output Point name
    Buffpath = os.path.join(os.getcwd(), "Buff")

    driver = ogr.GetDriverByName("ESRI Shapefile")


    # output_layer.CreateField(ogr.FieldDefn("line1_id", ogr.OFTInteger))
    # output_layer.CreateField(ogr.FieldDefn("line2_id", ogr.OFTInteger))
    outofrange = 0

    shp1 = driver.Open(shp1)
    input_layer1 = shp1.GetLayer()
    spatial_ref = input_layer1.GetSpatialRef()

    shp2 = driver.Open(shp2)
    input_layer2 = shp2.GetLayer()
    k = 0
    for feature1 in input_layer1:
        geom1 = feature1.GetGeometryRef()

        for feature2 in input_layer2:
            geom2 = feature2.GetGeometryRef()

            if geom1.Intersects(geom2):
                intersection = geom1.Intersection(geom2)
                if intersection is not None:
                    if intersection.ExportToWkt().split()[0] == "MULTIPOINT":
                        for p in intersection.ExportToWkt().split("(")[1][:-1].split(","):
                            Wkt = "POINT (" + p + ")"
                            BuffPointfile = os.path.join(Buffpath,
                                                         "BuffPoint" + str(int(Buffdistance)) + "_" + str(k) + ".shp")
                            output_source = driver.CreateDataSource(BuffPointfile)
                            output_layer = output_source.CreateLayer("intersection", spatial_ref, ogr.wkbPoint)
                            output_feature = ogr.Feature(output_layer.GetLayerDefn())
                            geom = ogr.CreateGeometryFromWkt(Wkt)
                            output_feature.SetGeometry(geom)
                            output_layer.CreateFeature(output_feature)
                            Searching_Points.append(BuffPointfile)
                            k += 1
                    else:
                        BuffPointfile = os.path.join(Buffpath, "BuffPoint" + str(int(Buffdistance)) + "_" + str(k) +".shp")
                        output_source = driver.CreateDataSource(BuffPointfile)
                        output_layer = output_source.CreateLayer("intersection", spatial_ref, ogr.wkbPoint)
                        output_feature = ogr.Feature(output_layer.GetLayerDefn())
                        geom = ogr.CreateGeometryFromWkt(intersection.ExportToWkt())
                        output_feature.SetGeometry(geom)
                        output_layer.CreateFeature(output_feature)
                        Searching_Points.append(BuffPointfile)
                        k += 1

    input_layer1 = input_layer2 = output_layer = geom1 = geom2 = None

    return Searching_Points

def CalcAreaVolume(Carray, Cmeta, PointDEM, drop, resolution):

    Copyarray = copy.deepcopy(Carray)
    Copyarray[Copyarray == np.nan] = -9999
    Copyarray[Copyarray < PointDEM] = np.nan
    Copyarray[Copyarray > PointDEM + drop] = np.nan
    Copyarray[Copyarray == -9999] = np.nan
    temparray = (np.ones([Cmeta['array_rows'], Cmeta['array_cols']]) * (PointDEM + drop)) - Copyarray
    temparray[Copyarray == np.nan] = np.nan
    Index = np.where(temparray > 0)
    Aream2 = abs(len(Index[0]) * (resolution ** 2))
    Volumem3 = np.nansum(temparray * (resolution ** 2))
    Areakm2 = abs(len(Index[0]) * (resolution ** 2)) / 1000000

    return Aream2, Areakm2, Volumem3, temparray

def Calc(Fname, DEM, wshp, pointnetshp, PointDEM, TargetHeight, resolution, delta):
    Err = 0
    ####
    Cutarray, CutMeta = warp(wshp, DEM)
    ####
    Datadic = {}

    tempDNS = FindDNStreamline(os.path.join(os.getcwd(), "temp\\" + Fname + "_tree.dat"))
    StreamDic = extractValue(os.path.join(os.getcwd(), "temp\\" + Fname + "_net.shp"))
    Maxdistance = 0
    Drop = 0

    for DNS in tempDNS:
        distance = 0
        drop = 0
        for linkno in DNS:
            distance += StreamDic[linkno][0]
            drop += StreamDic[linkno][1]
        if distance > Maxdistance:
            Maxdistance = distance
            Drop = drop
            Maxdistancelist = DNS

    Maxdistance = MakeLongStream(Maxdistancelist, pointnetshp, r".\temp.\LLLL.shp")

    Slope = Drop / (Maxdistance * 111.320 * delta * 1000)

    if Slope >= 0.005:
        Tc = 0.833 * ((Maxdistance * 111.320 * delta) / (Slope ** 0.6))
    else:
        Tc = 0.444 * ((Maxdistance * 111.320 * delta) / (Slope ** 0.515))

    drop = float(TargetHeight)

    Am2, Akm2, Volm3, TA = CalcAreaVolume(Cutarray, CutMeta, PointDEM, drop, resolution)

    if drop <= 0:
        fopen = open(r".\Result.\ErrorOccured.txt")
        fopen.write("Can't calculate DAM properties : Drop < 0")
        fopen.close()
        Err = 1
    else:
        while drop > 0:
            Am2, Akm2, Volm3, TTA = CalcAreaVolume(Cutarray, CutMeta, PointDEM, drop, resolution)

            Datadic[drop] = [Am2, Akm2, Volm3]
            drop -= 1

    return Datadic, Maxdistance, Slope, Tc, TA, CutMeta, Maxdistancelist, Err

def array2raster(newRasterfn, array, metadata):
    cols = array.shape[1]
    rows = array.shape[0]
    originX = metadata["geotransform"][0]
    originY = metadata["geotransform"][3]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, metadata["pixelWidth"], 0, originY, 0, metadata["pixelHeight"]))
    outRaster.SetProjection(metadata['projection'])
    outband = outRaster.GetRasterBand(1)
    if not metadata['noDataValue'] == None:
        outband.SetNoDataValue(metadata['noDataValue'])
    outband.WriteArray(array)
    outband.FlushCache()

def FindEndStreamline(treedat):
    EndStreamList = []
    with open(treedat) as Fopen:
        Flines = Fopen.read().splitlines()
        for k in range(0, len(Flines)):
            if int(Flines[k].split("\t")[5]) == -1 and int(Flines[k].split("\t")[6]) == -1:
                EndStreamList.append(int(Flines[k].split("\t")[1]))

    Streams = []
    for LN in EndStreamList:
        Streams.append(GetDNStreamNum(LN, treedat))

    return Streams

def GetMaxPotentialDrop(treedat, netdbf):
    streamlists = FindEndStreamline(treedat)
    with dbf.Table(netdbf) as CT:
        dbf.export(CT)
    dropdic = {}
    with open(netdbf[:-4] + ".csv", 'r') as fopen:
        flines = fopen.read().splitlines()
        for line in flines[1:]:
            datas = line.split(",")
            dropdic[int(datas[0].replace("\"", "", 2))] = float(datas[9].replace("\"", "", 2))
    MaxPoDrop = 0
    tempPoDrop = 0
    for streamlist in streamlists:
        for k in streamlist:
            tempPoDrop += dropdic[k]
        if tempPoDrop >= MaxPoDrop:
            MaxPoDrop = tempPoDrop
        tempPoDrop = 0

    return MaxPoDrop

def makeplot(diction, Maxdis, SL, Tc, TA, Cmeta, pointnetshp, MDL, resolution, delta, PointinfoDic, Pname):
    droplist = []
    Arealist = []
    Vollist = []
    keys = list(diction.keys())
    for key in keys:
        droplist.append(float(key))
        Arealist.append(diction[key][1])
        Vollist.append(diction[key][2])
    fig, axes = plt.subplots(1,2)
    axes[0].plot(droplist, Arealist, marker='o', color='black')
    axes[0].set_ylabel('$Area(km^{2}$)')
    axes[0].set_xlabel("Drop(m)")
    axes[0].grid(True)
    axes[1].plot(droplist, np.divide(Vollist, 1000000), marker='+', color='black')
    axes[1].set_ylabel('$Volume(10^{6} m^{3}$)')
    axes[1].set_xlabel("Drop(m)")
    axes[1].grid(True)
    fig.tight_layout()
    plt.savefig(".\Result.\AreaVol_"+ Pname + ".jpg")
    plt.clf()
    ##
    PointinfoDic[Pname]["Droplist"] = droplist
    PointinfoDic[Pname]["Arealist"] = Arealist # km2
    PointinfoDic[Pname]["Vollist"] = Vollist
    ##
    TL = MakeLongStream(MDL, pointnetshp, r".\temp.\LLL_" + Pname + ".shp")
    array2raster(r".\temp.\WHArea_"+ Pname + ".tif", TA, Cmeta)

    polygonize(r".\temp.\WHArea_"+ Pname + ".tif", r".\temp.\WHADArea_"+ Pname + ".shp")
    # dissolve(r".\temp.\WHADArea_"+ Pname + ".shp", r".\Result.\D_WHADArea_"+ Pname + ".shp") ##
    Length = intersects_L(r".\temp.\LLL_"+ Pname + ".shp", r".\temp.\WHADArea_"+ Pname + ".shp", r".\temp.\TLL_"+ Pname + ".shp")

    # Length = 1

    sorted_dic = sorted(diction.items(), key=operator.itemgetter(0))
    try:
        HperL = float(sorted_dic[-1][0]) / (Length * 111.320 * delta * 1000)
        VperA = float(sorted_dic[-1][1][2]) / float(sorted_dic[-1][1][0])
    except:
        HperL = 0
        VperA = 0

    PointinfoDic[Pname]["Length"] = (Maxdis * 111.320 * delta * 1000)
    PointinfoDic[Pname]["Slope"] = SL
    PointinfoDic[Pname]["Tc"] = Tc
    PointinfoDic[Pname]["H/L"] = HperL
    PointinfoDic[Pname]["V/A"] = VperA

def write_shapefile(poly, out_shp, spatial_ref):

    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(out_shp)
    layer = ds.CreateLayer('', spatial_ref, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
    defn = layer.GetLayerDefn()

    feat = ogr.Feature(defn)
    feat.SetField('id', 123)

    feat.SetGeometry(poly)

    layer.CreateFeature(feat)
    feat = geom = None

    ds = layer = feat = geom = None

def write_shapelinefile(pointlist, pointshp, out_shp):

    driver = ogr.GetDriverByName('Esri Shapefile')

    shp1 = driver.Open(pointshp)
    input_layer1 = shp1.GetLayer()
    spatial_ref = input_layer1.GetSpatialRef()

    ds = driver.CreateDataSource(out_shp)
    layer = ds.CreateLayer('', spatial_ref, ogr.wkbLineString)

    layer.CreateField(ogr.FieldDefn('id', ogr.OFTString))
    defn = layer.GetLayerDefn()

    feat = ogr.Feature(defn)
    feat.SetField('id', "Verticalline")
    line = ogr.Geometry(ogr.wkbLineString)
    line.AddPoint(float(pointlist[0][0]), float(pointlist[0][1]))
    line.AddPoint(float(pointlist[1][0]), float(pointlist[1][1]))

    feat.SetGeometry(line)

    layer.CreateFeature(feat)

    feat = geom = None

    ds = layer = feat = geom = None

def write_shapepointfile(pointdic, out_shp):

    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(out_shp)
    layer = ds.CreateLayer('', None, ogr.wkbMultiPoint)
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTString))
    defn = layer.GetLayerDefn()
    feat = ogr.Feature(defn)
    feat.SetField('id', "Overflowpoint")
    points = ogr.Geometry(ogr.wkbMultiPoint)
    keylist = list(pointdic.keys())
    for key in keylist:
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(float(key[0]), float(key[1]))
        points.AddGeometry(point)
    feat.SetGeometry(points)
    layer.CreateFeature(feat)
    feat = geom = None

    ds = layer = feat = geom = None

def GetDAMLine(X, Y, dis, dir, outshpname, pointshp):

    x = symbols('x')

    if dir == 7 or dir == 3:
        eqn = Eq(sqrt((X-x)**2), dis)
        Ans = solve(eqn, x)
        Point = [[Ans[0], Y], [Ans[1], Y]]


    if dir == 8 or dir == 4:
        b = Y - X
        eqn = Eq(sqrt((X-x)**2 + (Y-(x+b))**2), dis)
        Ans = solve(eqn, x)
        Point = [[Ans[0], Ans[0]+b], [Ans[1], Ans[1]+b]]

    if dir == 6 or dir == 2:
        b = Y + X
        eqn = Eq(sqrt((X - x) ** 2 + (Y - (-x + b)) ** 2), dis)
        Ans = solve(eqn, x)
        Point = [[Ans[0], (-Ans[0])+b], [Ans[1], (-Ans[1])+b]]

    if dir == 5 or dir == 1:
        eqn = Eq(sqrt((Y-x)**2), dis)
        Ans = solve(eqn, x)
        Point = [[X, Ans[0]], [X, Ans[1]]]

    write_shapelinefile(Point, pointshp, outshpname)

def FeatureToRaster(InputVector, RefImage, OutputImage):

    gdalformat = 'GTiff'
    datatype = gdal.GDT_Int16

    Image = gdal.Open(RefImage, gdal.GA_ReadOnly)


    Shapefile = ogr.Open(InputVector)
    Shapefile_layer = Shapefile.GetLayer()

    Output = gdal.GetDriverByName(gdalformat).Create(OutputImage, Image.RasterXSize, Image.RasterYSize, 1, datatype)
    Output.SetProjection(Image.GetProjectionRef())
    Output.SetGeoTransform(Image.GetGeoTransform())

    Band = Output.GetRasterBand(1)
    Band.SetNoDataValue(-9999)
    gdal.RasterizeLayer(Output, [1], Shapefile_layer)

    Band = None
    Output = None
    Image = None
    Shapefile = None

def pixelOffset2coord(raster, xOffset,yOffset):
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    coordX = originX+pixelWidth*xOffset
    coordY = originY+pixelHeight*yOffset
    return coordX, coordY

def AllPointList(DEM, TargetFeature):

    ELRasterOpen = gdal.Open(DEM)
    gt = ELRasterOpen.GetGeoTransform()
    ELRaster = ELRasterOpen.GetRasterBand(1)
    DEMdataarray = ELRaster.ReadAsArray()

    raster = gdal.Open(TargetFeature)
    Rasterband = raster.GetRasterBand(1)
    array = Rasterband.ReadAsArray()

    Pointlist = []

    points = np.where(array != -9999)
    for i in range(0, len(points[0])):
        Xcoord, Ycoord = pixelOffset2coord(raster, points[0][i], points[1][i])
        Pointlist.append(DEMdataarray[int((Xcoord - gt[0]) / gt[1]), int((Ycoord - gt[3]) / gt[5])])

    array = None
    DEMdataarray = None

    return(Pointlist)

def FindInflection(list, Resolution, VLD):
    pointlist = []

    if list[0] - list[1]< 0:
        Bval = "-"
    elif list[0] - list[1]> 0:
        Bval = "+"
    else:
        Bval = "="

    for k in range(0, len(list)-1):
        if list[k] - list[k+1] < 0:
            tempval = "-"
        elif list[k] - list[k+1] > 0:
            tempval = "+"
        else:
            tempval = "="
        if Bval != tempval:
            Bval = tempval
            pointlist.append(int(k))
    cen = np.asarray(pointlist) - int(len(list)/2)
    L = list[:int(len(list)/2)]
    R = list[int(len(list)/2):]
    temlistL = []
    temlistR = []
    for kk in cen:
        if kk <= 0:
            temlistL.append(list[kk + int(len(list)/2)])
        else:
            temlistR.append(list[kk + int(len(list)/2)])
    if len(temlistL) <= 1 and len(temlistR) <= 1:
        Suggestlinepo = min(list[0], list[-1])
    elif len(temlistL) <= 1:
        Suggestlinepo = max(temlistR)
    elif len(temlistR) <= 1:
        Suggestlinepo = max(temlistL)
    else:
        maxL = max(temlistL)
        maxR = max(temlistR)
        Suggestlinepo = min(maxL, maxR)

    if len(temlistL[-7:-1]) > 0 and len(temlistR[:6])>0:
        unline = min(min(temlistL[-7:-1]), min(temlistR[:6]))
    elif len(temlistL[-7:-1]) == 0 and len(temlistR[:6]) > 0:
        unline = min(temlistR[:6])
    elif len(temlistL[-7:-1]) > 0 and len(temlistR[:6]) == 0:
        unline = min(temlistL[-7:-1])
    else:
        min(list)

    leftside = [i for i, x in enumerate(L) if x >= Suggestlinepo]
    rightside = [i for i, x in enumerate(R) if x >= Suggestlinepo]

    if len(leftside) > 0 and len(rightside) > 0:
        leftdis = len(L) - max(leftside)
        rightdis = min(rightside)
        dis = leftdis + rightdis

    elif len(leftside) > 0 and len(rightside) == 0:
        dis = len(L) - max(leftside) + (VLD/2/Resolution)

    elif len(leftside) == 0 and len(rightside) > 0:
        dis = min(rightside) + (VLD/2/Resolution)
    else:
        dis = (VLD/Resolution)

    return pointlist, Suggestlinepo, unline, dis

def WarpAndGetData(inshp, array, meta, point, dis, TargetH, UnH):

    shp = ogr.Open(inshp)
    layer = shp.GetLayer()
    union = ogr.Geometry(ogr.wkbPolygon)
    gt = meta['geotransform']

    for feat in layer:
        geom = feat.GetGeometryRef()
        union = union.Union(geom)

    newuni = union.Buffer(0)

    Valuedic = {}

    BC = GenerateBuff(point, dis)

    Pointlist = Diff(BC, newuni)
    for Point in Pointlist:
        x = Point.split()[0]
        y = Point.split()[1]
        if x[0] == "(":
            x = x[1:]
        if y[-1] == ")":
            y = y[0:-1]
        if float(x) > 180 or float(x) < -180:
            continue
        if float(x) > 180 or float(x) < -180:
            continue
        mx, my = float(x), float(y)

        px = int((mx - gt[0]) / gt[1])
        py = int((my - gt[3]) / gt[5])

        Valuedic[float(x), float(y)] = float(array[py, px])

    Valkeys = list(Valuedic.keys())
    for dickey in Valkeys:
        if Valuedic[dickey] >= TargetH:
            del(Valuedic[dickey])
    return Valuedic

def MakeLongStream(MDL, pointstreamshp, outshp):

    driver = ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource(outshp)
    dslayer = ds.CreateLayer('', None, ogr.wkbMultiLineString)
    featureDefn = dslayer.GetLayerDefn()

    fieldNames = []

    shp1 = ogr.Open(pointstreamshp)
    layer1 = shp1.GetLayer()

    for i in range(layer1.GetLayerDefn().GetFieldCount()):
        fieldDefn = layer1.GetLayerDefn().GetFieldDefn(i)
        dslayer.CreateField(fieldDefn)
        fieldNames.append(fieldDefn.name)

    Length = 0
    for feature in layer1:
        if feature.GetField("LINKNO") in MDL:
            ingeom = feature.GetGeometryRef()
            Length += ingeom.Length()
            fieldVals = []
            for f in fieldNames:
                fieldVals.append(feature.GetField(f))
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(ingeom)

            for v, val in enumerate(fieldVals):
                outFeature.SetField(fieldNames[v], val)
            dslayer.CreateFeature(outFeature)

    return Length

def GenerateBuff(point, dis):
    Buffdistance = dis
    ds = ogr.Open(point)
    driver = ds.GetDriver()

    Buffcirclefile = os.path.join(temppath, "Buff.shp")
    if os.path.exists(Buffcirclefile):
        driver.DeleteDataSource(Buffcirclefile)
    outputbufferd = driver.CreateDataSource(Buffcirclefile)
    bufferlayer = outputbufferd.CreateLayer(Buffcirclefile, geom_type=ogr.wkbPolygon)
    featureDefn = bufferlayer.GetLayerDefn()

    fieldNames = []

    inlayer = ds.GetLayer()

    for i in range(inlayer.GetLayerDefn().GetFieldCount()):
        fieldDefn = inlayer.GetLayerDefn().GetFieldDefn(i)
        bufferlayer.CreateField(fieldDefn)
        fieldNames.append(fieldDefn.name)

    for feature in inlayer:
        ingeom = feature.GetGeometryRef()
        fieldVals = []
        for f in fieldNames: fieldVals.append(feature.GetField(f))
        outFeature = ogr.Feature(featureDefn)
        geomBuffer = ingeom.Buffer(Buffdistance)
        outFeature.SetGeometry(geomBuffer)

        for v, val in enumerate(fieldVals):
            outFeature.SetField(fieldNames[v], val)
        bufferlayer.CreateFeature(outFeature)
        outFeature = None

    return Buffcirclefile

def Diff(shp1, union):
    shp1 = ogr.Open(shp1)
    layer1 = shp1.GetLayer()

    union1 = ogr.Geometry(layer1.GetGeomType())
    for feat in layer1:
        geom = feat.GetGeometryRef()
        union1 = union1.Union(geom)

    union = union.GetBoundary()
    intersection = union.Difference(union1)

    Pointlist = intersection.ExportToWkt()[18:-2].split(",")

    feat = geom = None

    ds = layer = feat = geom = None

    return Pointlist

def DegreeToMeter(Rasterfile):
    metadata = {}
    dataset = gdal.Open(Rasterfile)
    metadata['array_rows'] = dataset.RasterYSize
    metadata['array_cols'] = dataset.RasterXSize
    metadata['geotransform'] = dataset.GetGeoTransform()
    metadata['projection'] = dataset.GetProjection()

    Lat = (metadata['geotransform'][3] + dataset.RasterYSize * metadata['geotransform'][5] + metadata['geotransform'][3])/2

    delta = np.cos(Lat * np.pi / 180)

    resolution = abs(metadata['geotransform'][1]) * 111.320 * delta * 1000

    metadata = dataset = Lat = None

    return resolution, delta

def Step1(Setdic, RiverDirection, Verlinelength, resolution, Suggestline, Unline, list2, GMPD, disk_30, point, pointshed, pointnetshp, delta, PointinfoDic, Pname=None):

    if Pname == None:
        Pname = "Start"

    PointinfoDic[Pname] = {}

    if RiverDirection in [2, 4, 6, 8]:
        dist = resolution * sqrt(2)
        disk_30 = float(disk_30 * dist)
        x2 = np.arange(int(-Verlinelength / 2), int(Verlinelength / 2), dist)
        if len(list2) > len(x2):
            x2 = x2.tolist()
            while len(list2) != len(x2):
                x2.append(x2[-1] + dist)
        elif len(list2) < len(x2):
            x2 = x2.tolist()
            while len(list2) != len(x2):
                x2 = x2[:-1]

    else:
        dist = resolution
        disk_30 = float(disk_30 * dist)
        x2 = np.arange(int(-Verlinelength / 2), int(Verlinelength / 2) + 1, dist)
        if len(list2) > len(x2):
            x2 = x2.tolist()
            while len(list2) != len(x2):
                x2.append(x2[-1] + dist)
        elif len(list2) < len(x2):
            x2 = x2.tolist()
            while len(list2) != len(x2):
                x2 = x2[:-1]

    disk = (disk_30/111.320) / delta / 1000

    texts = []
    TargetHeight = Suggestline
    Oarray, OMeta = raster2array(DEM)
    valdic = WarpAndGetData(pointshed, Oarray, OMeta, point, disk / 2, TargetHeight, Unline)

    try:
        Nooverflowline = min(valdic.values())
        if Nooverflowline <= Unline:
            Nooverflowline = TargetHeight

        L = list2[:int(len(list2) / 2)]
        R = list2[int(len(list2) / 2):]
        leftside = [i for i, x in enumerate(L) if x >= Nooverflowline]
        rightside = [i for i, x in enumerate(R) if x >= Nooverflowline]

        if len(leftside) > 0 and len(rightside) > 0:
            leftdis = len(L) - max(leftside)
            rightdis = min(rightside)
            dis = leftdis + rightdis

        elif len(leftside) > 0 and len(rightside) == 0:
            dis = len(L) - max(leftside) + (Verlinelength / 2 / resolution)

        elif len(leftside) == 0 and len(rightside) > 0:
            dis = min(rightside) + (Verlinelength / 2 / resolution)
        else:
            dis = (Verlinelength / resolution)

        disk = (float(dis * dist) / 111.320) / delta / 1000

    except:
        Nooverflowline = TargetHeight

    if Nooverflowline > GMPD + Unline:
        for po in [GMPD + Unline, Suggestline, Nooverflowline]:
            texts.append(plt.text(-4000, po, po))
        plt.axhline(y=Nooverflowline, color='r', linewidth=2)
        plt.axhline(y=Suggestline, color='b', linewidth=2)
        plt.axhline(y=GMPD + Unline, color='y', linewidth=2)
        plt.plot(x2, list2, color='grey')
        adjust_text(texts)
        plt.minorticks_on()
        plt.grid(True, which='minor', lw=0.6, ls="--", c=".90")
        plt.grid(True, which='major')
        plt.xlabel("Distance(m)")
        plt.ylabel("E.L.(m)")
        plt.legend(["E.L. : " + str(Nooverflowline) + " Drop : " + str(Nooverflowline - Unline),
                    "E.L. : " + str(Suggestline) + " Drop : " + str(Suggestline - Unline),
                    "E.L. : " + str(GMPD + Unline) + " Drop : " + str(GMPD)], loc='upper right')
        # Drop = Nooverflowline - Unline
        Drop = GMPD
        plt.savefig(".\Result.\Streamsection_PD_" + Pname + ".jpg")
        plt.clf()

    else:
        for po in [Suggestline, Nooverflowline]:
            texts.append(plt.text(-4000, po, po))
        plt.axhline(y=Nooverflowline, color='r', linewidth=2)
        plt.axhline(y=Suggestline, color='b', linewidth=2)
        plt.plot(x2, list2, color='grey')
        adjust_text(texts)
        plt.minorticks_on()
        plt.grid(True, which='minor', lw=0.6, ls="--", c=".90")
        plt.grid(True, which='major')
        plt.xlabel("Distance(m)")
        plt.ylabel("E.L.(m)")
        plt.legend(["E.L. : " + str(Nooverflowline) + " Drop : " + str(Nooverflowline - Unline),
                    "E.L. : " + str(Suggestline) + " Drop : " + str(Suggestline - Unline)], loc='upper right')
        Drop = Nooverflowline - Unline
        plt.savefig(".\Result.\Streamsection_" + Pname + ".jpg")
        plt.clf()

    dic, Maxdis, Sl, Tc, TA, Cmeta, MDL, Err = Calc(Fname, DEM, pointshed, pointnetshp, Unline, Drop, resolution, delta)

    if Err != 1:
        makeplot(dic, Maxdis, Sl, Tc, TA, Cmeta, pointnetshp, MDL, resolution, delta, PointinfoDic, Pname)

    PointinfoDic[Pname]["TargetDrop"] = Drop
    PointinfoDic[Pname]["Damwidth"] = disk * 111.320 * delta * 1000
    PointinfoDic[Pname]["PointHeight"] = Unline

    if Setdic['LAO'] == "1":
        LUdict = LUAREA(r".\temp.\WHADArea_"+ Pname + ".shp", r".\Supplementary.\MCD12Q1.tif")
        PointinfoDic[Pname]["LAO"] = LUdict
    if Setdic['SAO'] == "1":
        SAdict = SA(point, r".\temp.\WHADArea_"+ Pname + ".shp", r".\Supplementary.\Soil_Total.tif")
        PointinfoDic[Pname]["SAO"] = SAdict
    if Setdic['ASA'] == "1":
        Slope_D = GetValueFromTif( r".\Data\\" + Setdic['InFName'] + "_Slope_D.tif", point)
        PointinfoDic[Pname]["Slope_D"] = Slope_D[0]
    # If Option!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return PointinfoDic

def intersectsVicgrid(shp1, shp2, delta):

    rfopen = open(r".\Result.\Vic_grid_clip_Result.txt", 'w')
    sqDegreetosqKm = (111.320 * delta) ** 2
    shp1 = ogr.Open(shp1)
    layer1 = shp1.GetLayer()

    shp2 = ogr.Open(shp2)
    layer2 = shp2.GetLayer()

    fieldNames = []

    for i in range(layer1.GetLayerDefn().GetFieldCount()):
        fieldDefn = layer1.GetLayerDefn().GetFieldDefn(i)
        fieldNames.append(fieldDefn.name)

    union2 = ogr.Geometry(layer2.GetGeomType())

    for feat in layer2:
        geom = feat.GetGeometryRef()
        union2 = union2.Union(geom)

    Unigeom = ogr.ForceToMultiPolygon(union2)
    totarea = 0

    for feature1 in layer1:
        geom1 = feature1.GetGeometryRef()
        fieldVals = []
        for f in fieldNames:
            fieldVals.append(feature1.GetField(f))
        if geom1.Intersects(Unigeom):
            Intersection = Unigeom.Intersection(geom1)
            area = (Intersection.GetArea()) * sqDegreetosqKm
            rfopen.write("{0} : {1}, Area : {2}\n".format(fieldNames[2], fieldVals[2], area))
            totarea += area

    layer = feat = geom = None

def GetEPSGNUM(shp):
    driver = ogr.GetDriverByName('Esri Shapefile')
    shpfile = driver.Open(shp)
    layer = shpfile.GetLayer()
    return int(layer.GetSpatialRef().ExportToWkt().split("\"")[-2])

def ReprojectLayer(targetshp, refshp):
    indriver = ogr.GetDriverByName('ESRI Shapefile')

    inEPSG = GetEPSGNUM(targetshp)
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.ImportFromEPSG(inEPSG)

    outEPSG = GetEPSGNUM(refshp)
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.ImportFromEPSG(outEPSG)

    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    inDataSet = indriver.Open(targetshp)
    inLayer = inDataSet.GetLayer()

    outshp = r".\temp.\reprojectedGrid.shp"

    outdriver = ogr.GetDriverByName('ESRI Shapefile')

    if os.path.exists(outshp):
        outdriver.DeleteDataSource(outshp)
    outDataSet = outdriver.CreateDataSource(outshp)

    outLayer = outDataSet.CreateLayer("reproj_layer", geom_type=ogr.wkbMultiPolygon)


    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    outLayerDefn = outLayer.GetLayerDefn()

    inFeature = inLayer.GetNextFeature()
    while inFeature:

        geom = inFeature.GetGeometryRef()

        geom.Transform(coordTrans)

        outFeature = ogr.Feature(outLayerDefn)

        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))

        outLayer.CreateFeature(outFeature)

        outFeature = None
        inFeature = inLayer.GetNextFeature()

    file = open(r".\temp.\reprojectedGrid.prj", 'w')
    file.write(outSpatialRef.ExportToWkt())
    file.close()

    inDataSet = None
    outDataSet = None

    return outshp

def LUAREA(inshp, inraster):

    shp = ogr.Open(inshp)
    layer = shp.GetLayer()
    spatial_ref = layer.GetSpatialRef()
    union = ogr.Geometry(ogr.wkbPolygon)

    geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
    for feature in layer:
        geomcol.AddGeometry(feature.GetGeometryRef())
    convexhull = geomcol.ConvexHull()

    newuni = convexhull.Buffer(0)
    write_shapefile(newuni, ".\\temp\\WHADT.shp", spatial_ref)

    outraster = os.path.join(os.getcwd(), "temp\\temp.tif")

    OutTile = gdal.Warp(outraster, inraster,
                        format="GTiff",
                        multithread=True,
                        cutlineDSName= ".\\temp\\WHADT.shp",
                        cropToCutline=True)
    OutTile = None

    os.remove(".\\temp\\WHADT.shp")
    Cutarray, CutMeta = raster2array(outraster)
    os.remove(outraster)
    value, counts = np.unique(Cutarray[~np.isnan(Cutarray)], return_counts=True)
    value = value.tolist()
    counts = counts.tolist()
    tot = sum(counts)

    LUdict = {}
    for type in range(0, 14):
        try:
            percent = counts[value.index(type)] / tot
            LUdict[type] = percent*100

        except:
            LUdict[type] = 0

    Cutarray = CutMeta = None
    return LUdict

def SA(point, inshp, inraster):

    TBfile = ".\\temp\\Buff_T.shp"
    # 대상 지점 포인트 읽기
    ds = ogr.Open(point)
    driver = ds.GetDriver()
    inlayer = ds.GetLayer()
    spatial_ref = inlayer.GetSpatialRef()

    outputbufferd = driver.CreateDataSource(TBfile)
    bufferlayer = outputbufferd.CreateLayer('', spatial_ref, geom_type=ogr.wkbPolygon)

    for i in range(inlayer.GetLayerDefn().GetFieldCount()):
        fieldDefn = inlayer.GetLayerDefn().GetFieldDefn(i)
        bufferlayer.CreateField(fieldDefn)
    # write_shapefile(newuni, ".\\temp\\WHADT.shp", spatial_ref)
    for feature in inlayer:
        geom = feature.GetGeometryRef()
        buffer = geom.Buffer(1000 / 111.320 / delta / 1000)
        # write_shapefile(buffer, TBfile, spatial_ref)
        output_feature = ogr.Feature(bufferlayer.GetLayerDefn())
        output_feature.SetGeometry(buffer)
        for i in range(feature.GetFieldCount()):
            output_feature.SetField(i, feature.GetField(i))
        bufferlayer.CreateFeature(output_feature)

    ds = None
    outputbufferd = None

    outraster = os.path.join(os.getcwd(), "temp\\temp.tif")

    OutTile = gdal.Warp(outraster, inraster,
                        format="GTiff",
                        multithread=True,
                        cutlineDSName=TBfile,
                        cropToCutline=True)
    OutTile = None

    os.remove(TBfile)
    Cutarray, CutMeta = raster2array(outraster)
    os.remove(outraster)
    value, counts = np.unique(Cutarray[~np.isnan(Cutarray)], return_counts=True)
    value = value.tolist()
    counts = counts.tolist()
    tot = sum(counts)

    SAdict = {}
    for type in range(0, 9):
        try:
            percent = counts[value.index(type)] / tot
            SAdict[type] = percent * 100

        except:
            SAdict[type] = 0

    Cutarray = CutMeta = None

    return SAdict

def SettingFileLoad():
    print("Read Setting File..")
    SettingValues = {}
    SettingFiles = open(os.path.join(os.getcwd(), "Setting.txt"), "r")
    SettingData = SettingFiles.read().splitlines()
    SettingFiles.close()
    for ii in range(0, len(SettingData)):
        if len(SettingData[ii].split("|")) > 1:
            Value = SettingData[ii].split("|")[0].strip()
            Key = SettingData[ii].split("|")[1].split(":")[0].strip()
            SettingValues[Key] = Value
            print("{0}: {1}".format(Key, Value))

    return SettingValues
    # self.SettingValues = SettingValues
    # self.sdt = float(SettingValues["Optsdt"])  # simulation time (sec)
    # self.SimTimeDelta = datetime.timedelta(seconds=int(SettingValues["Optsdt"]))

def dissolve(input, output, multipoly=False, overwrite=False):
    ds = ogr.Open(input)
    lyr = ds.GetLayer()
    out_ds, out_lyr = createDS(output, ds.GetDriver().GetName(), lyr.GetGeomType(), lyr.GetSpatialRef(), overwrite)
    defn = out_lyr.GetLayerDefn()
    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    for feat in lyr:
        if feat.geometry():
            feat.geometry().CloseRings()
            wkt = feat.geometry().ExportToWkt()
            multi.AddGeometryDirectly(ogr.CreateGeometryFromWkt(wkt))
    union = multi.UnionCascaded()
    if multipoly is False:
        for geom in union:
            poly = ogr.CreateGeometryFromWkb(geom.ExportToWkb())
            feat = ogr.Feature(defn)
            feat.SetGeometry(poly)
            out_lyr.CreateFeature(feat)
    else:
        out_feat = ogr.Feature(defn)
        out_feat.SetGeometry(union)
        out_lyr.CreateFeature(out_feat)
        out_ds.Destroy()
    ds.Destroy()
    return True

def Buffer(netshp, point, Init_Buff_Dis, Buff_interval, Max_Buffer_Dis, delta):

    if not os.path.isdir(os.path.join(os.getcwd(), "Buff")):
        os.mkdir(os.path.join(os.getcwd(), "Buff"))
    Buffpath = os.path.join(os.getcwd(), "Buff")

    Buffdistance = Init_Buff_Dis /111.320 / delta / 1000

    Searching_Points = []

    while Buffdistance < Max_Buffer_Dis /111.320 / delta / 1000:
        Buffcirclefile = os.path.join(Buffpath, "BuffCircle" + str(int(Buffdistance * 111.320 * delta * 1000)) + ".shp")

        # featureDefn = bufferlayer.GetLayerDefn()

        fieldNames = []
        # 대상 지점 포인트 읽기
        ds = ogr.Open(point)
        driver = ds.GetDriver()
        inlayer = ds.GetLayer()

        outputbufferd = driver.CreateDataSource(Buffcirclefile)
        bufferlayer = outputbufferd.CreateLayer(Buffcirclefile, geom_type=ogr.wkbLineString)

        for i in range(inlayer.GetLayerDefn().GetFieldCount()):
            fieldDefn = inlayer.GetLayerDefn().GetFieldDefn(i)
            bufferlayer.CreateField(fieldDefn)

        for feature in inlayer:
            geom = feature.GetGeometryRef()
            buffer = geom.Buffer(Buffdistance)
            output_feature = ogr.Feature(bufferlayer.GetLayerDefn())
            output_feature.SetGeometry(buffer)
            for i in range(feature.GetFieldCount()):
                output_feature.SetField(i, feature.GetField(i))
            bufferlayer.CreateFeature(output_feature)

        ds = None
        outputbufferd = None

        # In Def. Intersects BuffPoints Generated
        intersects(Searching_Points, netshp, Buffcirclefile, Buffdistance * 111.320 * delta * 1000)
        Buffdistance += Buff_interval /111.320 / delta / 1000

    return Searching_Points
def ResultWrite(Dic, SetDic, type=None):

    if type == "Start":
        P_SSC_file = open(".\Result.\P_Start_SSC.csv", "w")
        P_SSC_file.write("Drop,Flood_Area,Volume\n")
        for k in range(0, len(Dic['Start']['Droplist'])):
            P_SSC_file.write("{0},{1},{2}\n".format(Dic['Start']['Droplist'][k],
                                                    Dic['Start']['Arealist'][k],
                                                    Dic['Start']['Vollist'][k]))
        P_SSC_file.close()

        Tresultsfile = open(".\Result.\PointInfo.csv", "w")
        Tresultsfile.write("Point_Name,EL,Longest_Liver_Length,Slope,Tc,H/L,V/A,DamHeight,DamWidth")
        if SetDic['ASA'] == "1":
            Tresultsfile.write(",Slope_Degree\n")
        else: Tresultsfile.write("\n")

        Tresultsfile.write("Start,{0},{1},{2},{3},{4},{5},{6},{7}".format(Dic['Start']['PointHeight'],
                                                                            Dic['Start']['Length'],
                                                                            Dic['Start']['Slope'],
                                                                            Dic['Start']['Tc'],
                                                                            Dic['Start']['H/L'],
                                                                            Dic['Start']['V/A'],
                                                                            Dic['Start']['TargetDrop'],
                                                                            Dic['Start']['Damwidth']))
        if SetDic['ASA'] == "1":
            Tresultsfile.write(",{0}\n".format(Dic['Start']['Slope_D']))
        else: Tresultsfile.write("\n")
        Tresultsfile.close()

        if SetDic['LAO'] == "1":
            LUresultsfile = open(".\Result.\PointLU.csv", "w")
            LUresultsfile.write("Point,1,2,3,4,5,6,7,8,9,10,11,12,13\n")
            LUresultsfile.write('Start,')
            for key in Dic['Start']['LAO'].keys():
                LUresultsfile.write("{0},".format(Dic['Start']['LAO'][key]))
            LUresultsfile.write("\n")
            LUresultsfile.close()
        if SetDic['SAO'] == "1":
            SAresultsfile = open(".\Result.\PointSA.csv", "w")
            SAresultsfile.write("Point,1,2,3,4,5,6,7,8,9\n")
            SAresultsfile.write("Start,")
            for key in Dic['Start']['SAO'].keys():
                SAresultsfile.write("{0},".format(Dic['Start']['SAO'][key]))
            SAresultsfile.write("\n")
            SAresultsfile.close()

    else:
        Tresultsfile = open(".\Result.\PointInfo.csv", "a")
        for key in list(Dic.keys()):

            Tresultsfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8}".format("P_"+str(key).split("Point")[1],
                                                                              Dic[key]['PointHeight'],
                                                                            Dic[key]['Length'],
                                                                            Dic[key]['Slope'],
                                                                            Dic[key]['Tc'],
                                                                            Dic[key]['H/L'],
                                                                            Dic[key]['V/A'],
                                                                            Dic[key]['TargetDrop'],
                                                                            Dic[key]['Damwidth']))
            if SetDic['ASA'] == "1":
                Tresultsfile.write(",{0}\n".format(Dic[key]['Slope_D']))
            else:
                Tresultsfile.write("\n")

            Tresultsfile.close()

            if SetDic['LAO'] == "1":
                LUresultsfile = open(".\Result.\PointLU.csv", "a")
                LUresultsfile.write(str(key).split("Point")[1] + ",")
                for key2 in Dic[key]['LAO'].keys():
                    LUresultsfile.write("{0},".format(Dic[key]['LAO'][key2]))
                LUresultsfile.write("\n")
                LUresultsfile.close()

            if SetDic['SAO'] == "1":
                SAresultsfile = open(".\Result.\PointSA.csv", "a")
                SAresultsfile.write(str(key).split("Point")[1] + ",")
                for key2 in Dic[key]['SAO'].keys():
                    SAresultsfile.write("{0},".format(Dic[key]['SAO'][key2]))
                SAresultsfile.write("\n")
                SAresultsfile.close()

            P_SSC_file = open(".\Result.\P_" + str(key).split("Point")[1] + "_SSC.csv", "w")
            P_SSC_file.write("Drop,Flood_Area,Volume\n")
            for k in range(0, len(Dic[key]['Droplist'])):
                P_SSC_file.write("{0},{1},{2}\n".format(Dic[key]['Droplist'][k],
                                                        Dic[key]['Arealist'][k],
                                                        Dic[key]['Vollist'][k]))
            P_SSC_file.close()





if __name__=="__main__":
    t1= time.time()

    if not os.path.isdir(os.path.join(os.getcwd(), "temp")):
        os.mkdir(os.path.join(os.getcwd(), "temp"))
    temppath = os.path.join(os.getcwd(), "temp")
    if not os.path.isdir(os.path.join(os.getcwd(), "Result")):
        os.mkdir(os.path.join(os.getcwd(), "Result"))
    Resultpath = os.path.join(os.getcwd(), "Result")

    SetDic = SettingFileLoad()

    print("Watershed Segmentation Start...")
    ## Point Result Diction Generation
    PointinfoDic = {}

    DEM = os.path.join(".\Data", SetDic["InFName"] + ".tif")
    # DEM2 = os.path.join(".\Data\DB30_LV5", Shedname2 + ".tif")
    AD8 = os.path.join(".\Data", SetDic["InFName"] + "_AreaD8.tif")
    D8 = os.path.join(".\Data", SetDic["InFName"] + "_D8.tif")
    point = os.path.join(".\Point", SetDic["InFName"] + ".shp")

    resolution, delta = DegreeToMeter(DEM)

    Fname = Shed(DEM, AD8, D8, "Start", point)
    pointshed = os.path.join(temppath, Fname + "_shed.shp")
    polygonize(os.path.join(temppath, Fname + "_w.tif"), pointshed)
    pointnetshp = os.path.join(temppath, Fname + "_net.shp")

    ######Start Point Analysis#########

    print("First point feasibility analysis start")
    X, Y, RiverDirection = FindMajorDirection(D8, point)

    outshpname = "Vertical_line.shp"
    VerlinelengthD = (float(SetDic["VLL"]) / 111.320) / delta / 1000

    GetDAMLine(X, Y, VerlinelengthD / 2 , RiverDirection, os.path.join(temppath, outshpname), point)

    FeatureToRaster(os.path.join(temppath, outshpname), DEM, os.path.join(temppath, "rasterized.tif"))

    Plist2 = AllPointList(DEM, os.path.join(temppath, "rasterized.tif"))

    inflactionlist_30, Suggestline_30, Unline_30, disk_30 = FindInflection(Plist2, resolution, float(SetDic["VLL"]))

    GMPD = GetMaxPotentialDrop(os.path.join(temppath, SetDic["InFName"] + "_" + "Start" + "_tree.dat"),
                               os.path.join(temppath, SetDic["InFName"] + "_" + "Start" + "_net.dbf"))

    PointinfoDic = Step1(SetDic, RiverDirection, float(SetDic["VLL"]), resolution, Suggestline_30, Unline_30, Plist2, GMPD, disk_30, point, pointshed, pointnetshp, delta, PointinfoDic)

    ResultWrite(PointinfoDic, SetDic, "Start")

    ######Buff Point Generation#########

    Searching_Points = Buffer(pointnetshp, point, int(SetDic["SID"]),
                              int(SetDic["SDI"]),
                              int(SetDic["SMD"]),
                              delta)
    print("First point feasibility analysis Done")
    ######Buff Point Analysis##########
    print("Start analyzing upstream points...")

    i = 1

    for Bpoint in Searching_Points:
        try:
            PointinfoDic = {}
            Pname = Bpoint.split("\\")[-1].split(".")[0]
            Fname = Shed(DEM, AD8, D8, Pname, Bpoint)
            pointshed = os.path.join(temppath, Fname + "_shed.shp")
            polygonize(os.path.join(temppath, Fname + "_w.tif"), pointshed)
            pointnetshp = os.path.join(temppath, Fname + "_net.shp")
            X, Y, RiverDirection = FindMajorDirection(D8, Bpoint)
            GetDAMLine(X, Y, VerlinelengthD / 2, RiverDirection, os.path.join(temppath, outshpname), Bpoint)
            FeatureToRaster(os.path.join(temppath, outshpname), DEM, os.path.join(temppath, "rasterized.tif"))
            Plist2 = AllPointList(DEM, os.path.join(temppath, "rasterized.tif"))

            inflactionlist_30, Suggestline_30, Unline_30, disk_30 = FindInflection(Plist2, resolution, float(SetDic["VLL"]))
            GMPD = GetMaxPotentialDrop(os.path.join(temppath, SetDic["InFName"] + "_" + Pname + "_tree.dat"),
                                       os.path.join(temppath, SetDic["InFName"] + "_" + Pname + "_net.dbf"))
            PointinfoDic = Step1(SetDic, RiverDirection, float(SetDic["VLL"]), resolution, Suggestline_30, Unline_30, Plist2, GMPD,
                                 disk_30, Bpoint, pointshed, pointnetshp, delta, PointinfoDic, Pname)
            print("{0}/{1} Done".format(i, len(Searching_Points)))
            ResultWrite(PointinfoDic, SetDic)
        except:
            print("{0}/{1} Please Check The Point".format(i, len(Searching_Points)))

        i += 1

    print("Total Time Consumption: ", time.time() - t1)
