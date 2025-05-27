

// GDALtest.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include "include/gdal.h"
#include "include/gdal_priv.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "tinystr.h"
#include "tinyxml.h"
#include <cstdlib> // 包含atoi函数的头文件
using namespace std;
//控制点信息
struct controlPoint
{
    //获取时间
    string time;
    //x坐标
    double x;
    //y坐标
    double y;
    //z坐标
    double z;
    //x速度
    double vx;
    //y速度
    double vy;
    //z速度
    double vz;
};
//影像数据
struct metaData
{
    //影像行数
    int rows;
    //影像列数
    int cols;
    //轨道控制点数
    int controlPointsNum;
    //控制点信息
    vector<controlPoint> controlPoints;
    //距离视数
    int rangeLooks;
    //距离向空间分辨率
    double rowspasing;
    //方位向空间分辨率
    double azimuthRes;
    //载波中心频率
    double centerFrequency;
    //PRF
    double prf;
    //距离向采样频率
    double rsf;
    //成像最近距离
    double nearRange;
    //升降轨
    double orbitDirection;
    //侧视方向
    double lookDirection;
    //成像开始时间
    string startTime; 
};
//读入影像数据，输出影像长宽以及格式信息
GDALDataset* Read_image(GDALDataset* pDatasetRead, char* filename);
metaData Read_xml(char* filename);
//用于打印元数据的函数
void printMetaData(const metaData& md);
void printControlPoint(const controlPoint& cp);



int main()
{
    GDALAllRegister(); //注册GDAL库中的所有驱动程序
    GDALDataset* pDatasetRead{};
    GDALDataset* pDatasetSave;
    GDALDriver* pDriver;
    //char filename[200] = "data/imagery_HH.tif";
    char filename[200] = "data/t-intesnsity.tif";
    char xmlfilename[200] = "data/TSX1_SAR__SSC______SM_S_SRA_20121208T162059_20121208T162107.xml";

    pDatasetRead = Read_image(pDatasetRead, filename);
    //读取影像数据
    metaData data = Read_xml(xmlfilename);

    printMetaData(data);
   
    //释放内存和关闭数据集
    GDALClose(pDatasetRead);

}


GDALDataset* Read_image(GDALDataset* pDatasetRead, char* filename) {
    pDatasetRead = (GDALDataset*)GDALOpen(filename, GA_ReadOnly);
    if (pDatasetRead == NULL) {
        printf("Failed to open file: %s\n", filename);
        return pDatasetRead;
    }
    int lWidth = pDatasetRead->GetRasterXSize();
    int lHeight = pDatasetRead->GetRasterYSize();
    int nBands = pDatasetRead->GetRasterCount();
    GDALDataType datatype = pDatasetRead->GetRasterBand(1)->GetRasterDataType();
    cout << "width: " << lWidth << endl;
    cout << "height: " << lHeight << endl;
    cout << "Bands: " << nBands << endl;
    cout << "DataType: " << GDALGetDataTypeName(datatype) << endl;
    return pDatasetRead;
}
metaData Read_xml(char* filename) {
    metaData data;
    TiXmlDocument doc(filename);  // 读取XML文件
    bool loadOkay = doc.LoadFile();  // 判断是否载入成功
    if (!loadOkay) {
        cout << "Could not load XML file." << endl;
        return data;
    }
    // 获取根节点
    TiXmlElement* root = doc.FirstChildElement("level1Product");
    // 遍历节点路径
    //productInfo中获取影像行列数\成像最近距离、升降轨、侧视方向、成像开始时间
    TiXmlElement* productInfo = root->FirstChildElement("productInfo");
    TiXmlElement* imageDataInfo = productInfo->FirstChildElement("imageDataInfo");
    TiXmlElement* imageRaster = imageDataInfo->FirstChildElement("imageRaster");
    TiXmlElement* numberOfColumns = imageRaster->FirstChildElement("numberOfColumns");
    const char* numberofColumns = numberOfColumns->GetText();
    data.cols = stoi(numberofColumns);  // 获取列数
    TiXmlElement* numberOfRows = imageRaster->FirstChildElement("numberOfRows");
    const char* numberofRows = numberOfRows->GetText();
    data.rows = stoi(numberofRows);  // 获取行数
    TiXmlElement* sceneInfo = productInfo->FirstChildElement("sceneInfo");
    TiXmlElement* rangeTime = sceneInfo->FirstChildElement("rangeTime");
    TiXmlElement* firstPixel = rangeTime->FirstChildElement("firstPixel");
    const char* firstPixelValue = firstPixel->GetText();
    data.nearRange = atoi(firstPixelValue);  // 获取成像最近距离
    TiXmlElement* missionInfo = productInfo->FirstChildElement("missionInfo");
    TiXmlElement* orbitDirection = missionInfo->FirstChildElement("orbitDirection");
    const char* orbitDirectionValue = orbitDirection->GetText();
    if (strcmp(orbitDirectionValue, "ASCENDING") == 0) {
		data.orbitDirection = 1.0;  // 升轨
	} else if (strcmp(orbitDirectionValue, "DESCENDING") == 0) {
		data.orbitDirection = -1.0;  // 降轨
	} else {
		std::cerr << "Unknown orbit direction: " << orbitDirectionValue << std::endl;
		return data;
	}
    TiXmlElement* acquisitionInfo = productInfo->FirstChildElement("acquisitionInfo");
    TiXmlElement* lookDirection = acquisitionInfo->FirstChildElement("lookDirection");
    const char* lookDirectionValue = lookDirection->GetText();
    if (strcmp(orbitDirectionValue, "RIGHT") == 0) {
        data.lookDirection = 1.0;  // 右侧视
    }
    else if (strcmp(orbitDirectionValue, "DESCENDING") == 0) {
        data.lookDirection = -1.0;  // 左侧视
    }
    else {
        std::cerr << "Unknown orbit direction: " << orbitDirectionValue << std::endl;
        return data;
    }
    //读控制点信息
    TiXmlElement* platform = root->FirstChildElement("platform");
    TiXmlElement* orbit = platform->FirstChildElement("orbit");
    TiXmlElement* orbitHeader = orbit->FirstChildElement("orbitHeader");
    TiXmlElement* numStateVectors = orbitHeader->FirstChildElement("numStateVectors");
    // 获取节点的文本内容
    const char* numStateVectorsValue = numStateVectors->GetText();
    data.controlPointsNum = atoi(numStateVectorsValue);
    vector<controlPoint> cps;
    // 遍历所有 stateVec 节点
    for (TiXmlElement* stateVec = orbit->FirstChildElement("stateVec"); stateVec != nullptr; stateVec = stateVec->NextSiblingElement("stateVec")) {
        controlPoint cp;
        // 读取子节点
        cp.time = stateVec->FirstChildElement("timeUTC")->GetText();
        const char* posX = stateVec->FirstChildElement("posX")->GetText();
        cp.x = atof(posX);
        const char* posY = stateVec->FirstChildElement("posY")->GetText();
        cp.y = atof(posY);
        const char* posZ = stateVec->FirstChildElement("posZ")->GetText();
        cp.z = atof(posZ);
        const char* velX = stateVec->FirstChildElement("velX")->GetText();
        cp.vx = atof(velX);
        const char* velY = stateVec->FirstChildElement("velY")->GetText();
        cp.vy = atof(velY);
        const char* velZ = stateVec->FirstChildElement("velZ")->GetText();
        cp.vz = atof(velZ);

        cps.push_back(cp);
    }
    data.controlPoints = cps;
    TiXmlElement* processing = root->FirstChildElement("processing");
    TiXmlElement* processingParameter = processing->FirstChildElement("processingParameter");
    TiXmlElement* rangeLooks = processingParameter->FirstChildElement("rangeLooks");
    const char* rangeLooksValue = rangeLooks->GetText();
    data.rangeLooks = atoi(rangeLooksValue);//距离向视数
    TiXmlElement* rowSpacing = imageRaster->FirstChildElement("rowSpacing");
    const char* rowSpacingValue = rowSpacing->GetText();
    data.rowspasing = atof(rowSpacingValue); // 距离向空间分辨率
    TiXmlElement* calibration = root->FirstChildElement("calibration");
    TiXmlElement* nominalGeometricPerformance = calibration->FirstChildElement("nominalGeometricPerformance");
    TiXmlElement* azimuthRes = nominalGeometricPerformance->FirstChildElement("azimuthRes");
    const char* azimuthResValue = azimuthRes->GetText();
    data.azimuthRes = atof(azimuthResValue); // 方位向空间分辨率
    TiXmlElement* instrument = root->FirstChildElement("instrument");
    TiXmlElement* radarParameters = instrument->FirstChildElement("radarParameters");
    TiXmlElement* centerFrequency = radarParameters->FirstChildElement("centerFrequency");
    const char* centerFrequencyValue = centerFrequency->GetText();
    data.centerFrequency = atof(centerFrequencyValue); // 载波中心频率
    TiXmlElement* productSpecific = root->FirstChildElement("productSpecific");
    TiXmlElement* complexImageInfo = productSpecific->FirstChildElement("complexImageInfo");
    TiXmlElement* commonPRF = complexImageInfo->FirstChildElement("commonPRF");
    const char* commonPRFValue = commonPRF->GetText();
    data.prf = atof(commonPRFValue); // prf
    TiXmlElement* commonRSF = complexImageInfo->FirstChildElement("commonRSF");
    const char* commonRSFValue = commonRSF->GetText();
    data.rsf = atof(commonRSFValue); // rsf
    return data;
}

void printMetaData(const metaData& md) {
    cout << "Rows: " << md.rows << endl;
    cout << "Cols: " << md.cols << endl;
    cout << "Control Points Num: " << md.controlPointsNum << endl;
    cout << "Range Looks: " << md.rangeLooks << endl;
    cout << "Row Spacing: " << md.rowspasing << endl;
    cout << "Azimuth Resolution: " << md.azimuthRes << endl;
    cout << "Center Frequency: " << md.centerFrequency << endl;
    cout << "PRF: " << md.prf << endl;
    cout << "RSF: " << md.rsf << endl;
    cout << "Near Range: " << md.nearRange << endl;
    cout << "Orbit Direction: " << md.orbitDirection << endl;
    cout << "Look Direction: " << md.lookDirection << endl;
    cout << "Start Time: " << md.startTime << endl;

    cout << "Control Points:" << endl;
    for (const auto& cp : md.controlPoints) {
        printControlPoint(cp);
        cout << "------------------------" << endl;
    }
}
void printControlPoint(const controlPoint& cp) {
    cout << "Time: " << cp.time << endl;
    cout << "X: " << cp.x << endl;
    cout << "Y: " << cp.y << endl;
    cout << "Z: " << cp.z << endl;
    cout << "Vx: " << cp.vx << endl;
    cout << "Vy: " << cp.vy << endl;
    cout << "Vz: " << cp.vz << endl;
}