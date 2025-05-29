

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
#include "matrix.h"
using namespace std;
//控制点信息
struct controlPoint
{
    //获取时间
    int time;
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
    int startTime; 
    // 影像角点坐标
    vector<double> Lat; 
    vector<double> Lon; 
};
//读入影像数据，输出影像长宽以及格式信息
GDALDataset* Read_image(GDALDataset* pDatasetRead, char* filename);
metaData Read_xml(char* filename);
//用于打印元数据的函数
void printMetaData(const metaData& md);
void printControlPoint(const controlPoint& cp);
//轨道拟合
Matrix orbitFit(metaData data);
//求瞬时速度和加速度
Matrix GetSpeed(int t, Matrix A);
Matrix GetAddSpeed(int t, Matrix A);
//坐标转换
vector<Matrix> toXYZ(GDALDataset* pDatasetDEM, double dem_pixel, double lat0, double lon0);
//几何校正
void GeoCorrection(GDALDataset* pDatasetReadImg, vector<Matrix> dem);


int main()
{
    GDALAllRegister(); //注册GDAL库中的所有驱动程序
    GDALDataset* pDatasetReadImg{};
    GDALDataset* pDatasetDEM{};
    GDALDataset* pDatasetSave;
    GDALDriver* pDriver;
    //char filename[200] = "data/imagery_HH.tif";
    char Imgfilename[200] = "data/t-intensity.tif";
    char DEMfilename[200] = "data/dem_r.tif";
    char xmlfilename[200] = "data/TSX1_SAR__SSC______SM_S_SRA_20121208T162059_20121208T162107.xml";
    double dem_pixel = 4.63e-5;//加密后的像元大小（单位为度）
    //读入影像
    pDatasetReadImg = Read_image(pDatasetReadImg, Imgfilename);
    //读取影像数据
    metaData data = Read_xml(xmlfilename);
    //输出元数据
    printMetaData(data);
    //轨道拟合
    Matrix A = orbitFit(data);
    //求瞬间速度和加速度
    Matrix V = GetSpeed(100, A); // 计算t=0时刻的速度
    Matrix V_add = GetAddSpeed(100, A); // 计算t=0时刻的加速度
    //将DEM转为地心坐标
    pDatasetDEM = Read_image(pDatasetDEM, DEMfilename);//读入DEM数据
    double lat0 = 19.800139; // DEM左上角纬度:来自裁剪参数
    double lon0 = -155.700139; // DEM左上角经度
    vector<Matrix> DEMxyz = toXYZ(pDatasetDEM, dem_pixel, lat0, lon0); // 转换为地心坐标系

    //释放内存和关闭数据集
    GDALClose(pDatasetReadImg);

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
    cout << "--------------------------" << endl;
    cout << "成功打开文件：" << endl;
    cout << "File: " << filename << endl;
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
        const char* time = stateVec->FirstChildElement("timeGPS")->GetText();
        cp.time = atoi(time);  // 获取时间
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

    // 遍历所有 sceneCornerCoord 节点
    for (TiXmlElement* corner = sceneInfo->FirstChildElement("sceneCornerCoord");
        corner != nullptr;
        corner = corner->NextSiblingElement("sceneCornerCoord")) {


        // 提取纬度和经度
        double lat = atof(corner->FirstChildElement("lat")->GetText());
        data.Lat.push_back(lat); // 存储纬度
        double lon = atof(corner->FirstChildElement("lon")->GetText());
        data.Lon.push_back(lon); // 存储经度

    }

    return data;
}
void printMetaData(const metaData& meta) {
    cout << "Rows: " << meta.rows << endl;
    cout << "Cols: " << meta.cols << endl;
    cout << "Control Points Num: " << meta.controlPointsNum << endl;
    cout << "Range Looks: " << meta.rangeLooks << endl;
    cout << "Row Spacing: " << meta.rowspasing << endl;
    cout << "Azimuth Resolution: " << meta.azimuthRes << endl;
    cout << "Center Frequency: " << meta.centerFrequency << endl;
    cout << "PRF: " << meta.prf << endl;
    cout << "RSF: " << meta.rsf << endl;
    cout << "Near Range: " << meta.nearRange << endl;
    cout << "Orbit Direction: " << meta.orbitDirection << endl;
    cout << "Look Direction: " << meta.lookDirection << endl;
    cout << "Start Time: " << meta.startTime << endl;

    cout << "Control Points:" << endl;
    for (const auto& cp : meta.controlPoints) {
        printControlPoint(cp);
        cout << "------------------------" << endl;
    }
    cout << "sceneCornerCoord:" << endl;
    for (int i = 0; i < 4;i++) {
        cout << meta.Lat[i]<<"\t"<<meta.Lon[i] << endl;
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
Matrix orbitFit(metaData data) {
    Matrix T(data.controlPointsNum, 4); // T 矩阵
    Matrix A(4, 3); // 待求系数矩阵4行对应4个系数，3列对应x、y、z坐标
    Matrix Ax(4, 1); // Ax 矩阵
    Matrix Ay(4, 1); // Ay 矩阵
    Matrix Az(4, 1); // Az 矩阵
    for (int i = 0; i < data.controlPointsNum; i++) {
        cout << data.controlPoints[i].time - data.controlPoints[0].time << endl; // 输出时间差
        T.set(i, 0, 1.0); // 第一列为1
        T.set(i, 1, data.controlPoints[i].time - data.controlPoints[0].time); // 第二列为时间
        T.set(i, 2, (data.controlPoints[i].time - data.controlPoints[0].time) * (data.controlPoints[i].time - data.controlPoints[0].time)); // 第三列为x坐标
        T.set(i, 3, (data.controlPoints[i].time - data.controlPoints[0].time) * (data.controlPoints[i].time - data.controlPoints[0].time) * (data.controlPoints[i].time - data.controlPoints[0].time)); // 第四列为y坐标
    }
    cout << "T矩阵" << endl;
    T.print();
    Matrix X(data.controlPointsNum, 1); // X矩阵
    for (int i = 0; i < data.controlPointsNum; i++) {
        X.set(i, 0, data.controlPoints[i].x); // 设置X矩阵的值为控制点的x坐标
    }
    Ax = (T.t() * T).inv() * T.t() * X;
    cout << "A-x" << endl;
    Ax.print(); 
    Matrix L(data.controlPointsNum, 1);//残差值
    L = T * Ax - X; 
    cout << "残差-x" << endl;
    L.print(); 
    for (int i = 0; i < 4; i++) {
        A.set(i, 0, Ax.get(i, 0)); // 将Ax矩阵的值赋给A矩阵
    }
    for (int i = 0; i < data.controlPointsNum; i++) {
		X.set(i, 0, data.controlPoints[i].y); // 设置X矩阵的值为控制点的y坐标
	}
    Ay = (T.t() * T).inv() * T.t() * X;
    cout << "A-y" << endl;
    Ay.print();
    L = T * Ay - X;
    cout << "残差-y" << endl;
    L.print();
    for (int i = 0; i < 4; i++) {
		A.set(i, 1, Ay.get(i, 0)); // 将Ay矩阵的值赋给A矩阵
	}
    for (int i = 0; i < data.controlPointsNum; i++) {
        X.set(i, 0, data.controlPoints[i].z); // 设置X矩阵的值为控制点的z坐标
        }
    Az = (T.t() * T).inv() * T.t() * X;
    cout << "A-z" << endl;
    Az.print();
    L = T * Az - X;
    cout << "残差-z" << endl;
    L.print();
    for (int i = 0; i < 4; i++) {
		A.set(i, 2, Az.get(i, 0)); // 将Az矩阵的值赋给A矩阵
	}
    cout << "A矩阵" << endl;
    A.print();

    return A;
}

Matrix GetSpeed(int t, Matrix A)
{
    Matrix V(3, 1); // 速度矩阵代表xyz方向的速度
    for(int i=0;i<3;i++) // 对于每个方向的速度
	{
		V.set(i, 0, A.get(1, i) + 2 * A.get(2, i) * t + 3 * A.get(3, i) * t * t); // 速度公式
	}
    cout << "t时刻速度矩阵" << endl;
    V.print();
    cout<<"速度合成结果：" << sqrt(V.get(0, 0) * V.get(0, 0) + V.get(1, 0) * V.get(1, 0) + V.get(2, 0) * V.get(2, 0)) << endl; // 输出速度模
    return V;
}
Matrix GetAddSpeed(int t, Matrix A)
{
    Matrix V(3, 1); // 速度矩阵代表xyz方向的速度
    for (int i = 0; i < 3; i++) // 对于每个方向的速度
    {
        V.set(i, 0, 2 * A.get(2, i) + 6 * A.get(3, i) * t); // 速度公式
    }
    cout << "t时刻加速度矩阵" << endl;
    V.print();
    cout << "加速度合成结果：" << sqrt(V.get(0, 0) * V.get(0, 0) + V.get(1, 0) * V.get(1, 0) + V.get(2, 0) * V.get(2, 0)) << endl; // 输出速度模
    return V;
}
vector<Matrix> toXYZ(GDALDataset* pDatasetDEM, double dem_pixel, double lat0, double lon0) {
    vector<Matrix> result;//存储顺序为X、Y、Z
    Matrix X(pDatasetDEM->GetRasterYSize(), pDatasetDEM->GetRasterXSize()); // X坐标矩阵
    Matrix Y(pDatasetDEM->GetRasterYSize(), pDatasetDEM->GetRasterXSize()); // Y坐标矩阵
    Matrix Z(pDatasetDEM->GetRasterYSize(), pDatasetDEM->GetRasterXSize()); // Z坐标矩阵
    //DEM行列数读取
    int rows = pDatasetDEM->GetRasterYSize();
    int cols = pDatasetDEM->GetRasterXSize();
    //逐像素处理
    double lat = 0, lon = 0;
    for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
            //初始化
            double x = 0;
            double y = 0;
            double z = 0;
			//获取DEM像素值
			pDatasetDEM->GetRasterBand(1)->RasterIO(GF_Read, j, i, 1, 1, &z, 1, 1, GDT_Float64, 0, 0);
			//计算经纬度
			lat = lat0 + (i * dem_pixel); // 假设lat0为左上角纬度
			lon = lon0 + (j * dem_pixel); // 假设lon0为左上角经度
            // 椭球信息
            double a = 6378137.0; // 长半轴
            double b = 6356752.3142; // 短半轴
            double e2 = (a * a - b * b) / (a * a); // 第一偏心率平方
            // 将经纬度转换为弧度
            double lat_rad = lat * M_PI / 180.0; // 纬度转换为弧度
            double lon_rad = lon * M_PI / 180.0; // 经度转换为弧度
            // 计算地心坐标系的X、Y、Z坐标
            double N = a / sqrt(1 - e2 * sin(lat_rad) * sin(lat_rad)); // 卯酉圈半径
            x = (N + z) * cos(lat_rad) * cos(lon_rad); // 经度作为X坐标
            y = (N + z) * cos(lat_rad) * sin(lon_rad); // 纬度作为Y坐标
            z = ((1 - e2) * N + z) * sin(lat_rad); // 高程作为Z坐标
			//将经纬度转换为地心坐标系
			X.set(i, j, x); // 经度作为X坐标
			Y.set(i, j, y); // 纬度作为Y坐标
			Z.set(i, j, z); // 高程作为Z坐标
		}
	}
	result.push_back(X);
	result.push_back(Y);
	result.push_back(Z);
	return result; // 返回地心坐标系的X、Y、Z矩阵
}