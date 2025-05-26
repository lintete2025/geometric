

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
using namespace std;
//控制点信息
struct controlPoint
{
    //获取时间
    char time[20];
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
    double startTime; 
};
//读入影像数据，输出影像长宽以及格式信息
GDALDataset* Read_image(GDALDataset* pDatasetRead, char* filename);


int main()
{
    GDALAllRegister(); //注册GDAL库中的所有驱动程序
    GDALDataset* pDatasetRead{};
    GDALDataset* pDatasetSave;
    GDALDriver* pDriver;
    //char filename[200] = "data/imagery_HH.tif";
    char filename[200] = "data/t-intesnsity.tif";

    pDatasetRead = Read_image(pDatasetRead, filename);
    //读取影像数据
    TiXmlDocument doc("data/TSX1_SAR__SSC______SM_S_SRA_20121208T162059_20121208T162107.xml");  // 读取XML文件
    bool loadOkay = doc.LoadFile();  // 判断是否载入成功
    if (!loadOkay)
        cout << "Could not load XML file." << endl;
    //TiXmlElement* root = doc.FirstChildElement(); // 获取根节点
    // 获取根节点
    TiXmlElement* root = doc.FirstChildElement("level1Product");
    // 遍历节点路径
    TiXmlElement* platform = root->FirstChildElement("platform");
    TiXmlElement* orbit = platform->FirstChildElement("orbit");
    TiXmlElement* orbitHeader = orbit->FirstChildElement("orbitHeader");
    TiXmlElement* numStateVectors = orbitHeader->FirstChildElement("numStateVectors");
    // 获取节点的文本内容
    const char* numStateVectorsValue = numStateVectors->GetText();
    if (!numStateVectorsValue) {
        std::cerr << "Failed to get text from 'numStateVectors'." << std::endl;
        return -1;
    }

    // 输出结果
    std::cout << "轨道控制点数: " << numStateVectorsValue << std::endl;



   
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
