
/***********************************************
 * File: main.cpp
 *
 * Author: LYG
 * Date: °ËÔÂ 2018
 *
 * Purpose:
 *
 * 
 **********************************************/
#include <stdio.h>
#include <windows.h>
#include "Interface.h"
#include "FreeImage.h"
#include "Math/Math.h"
#include "Math/Quaternion.h"
#include "Math/Vector3.h"
#pragma comment(lib, "PanoramaDLL.lib")
#pragma comment(lib, "FreeImaged.lib")

#define SAFE_DELETE(p)       { if(p) { delete (p);     (p)=nullptr; } }
#define SAFE_DELETE_ARRAY(p) { if(p) { delete[] (p);   (p)=nullptr; } }
#define SAFE_RELEASE(p)      { if(p) { (p)->Release(); (p)=nullptr; } }

namespace ApplicationPath
{
	std::string CurrentPath;
	std::string GetApplicationPath()
	{
		if (CurrentPath.size() == 0)
		{
			char Temp[256];
			memset(Temp, 0, 256);
			GetModuleFileNameA(NULL, Temp, 256);
			int nLen = strlen(Temp);
			while (nLen > 0)
			{
				if (Temp[nLen] == '\\' || Temp[nLen] == '/')
				{
					break;
				}
				Temp[nLen--] = '\0';
			}
			CurrentPath = Temp;
		}
		return CurrentPath;
	}
}

void StartBuildPanorama()
{
	const int PanoramaWidth = 4096;
	const int PanoramaHeight = 2048;
	const int Width = 1920;
	const int Height = 1080;
	const int ImageCount = 7;
	std::string CurrentPath;
	
	float Fovy = 90;
	float RotateDegree = 72;
	float Near = 1.0f;
	float Far = 10000.0f;
	SetPanoramaRendererType(RendererD3D11);
	SetPanoramaWidthAndHeight(PanoramaWidth, PanoramaHeight);
	SetInputImageWidthAndHeight(Width, Height);
	SetInputImageCount(ImageCount);
	SetInputImageProjectionParameters(Fovy, Near, Far);

	float RotateRadian = Math::DegreesToRadians(360.0f / 5.0f);
	for (int k = 0; k < ImageCount; k++)
	{
		Quaternion q = Quaternion::IDENTITY;
		if (k == 5)
		{
			// left picture
			q.FromAngleAxis(Radian(Degree(90)), Vector3::UNIT_Y);
		}
		else if (k == 6)
		{
			// right picture
			q.FromAngleAxis(Radian(Degree(-90)), Vector3::UNIT_Y);
		}
		else if (k == 0)
		{
			// nothing to rotate
		}
		else
		{
			q.FromAngleAxis(Radian(-RotateRadian * k), Vector3::UNIT_X);
		}
		q.normalise();
		float ImageRotation[4] = {q.x, q.y, q.z, q.w};

		SetInputImageRotation(k, ImageRotation);
	}

	Initialise();
	CurrentPath = ApplicationPath::GetApplicationPath();
	// Initialise Input Texture Data
	unsigned char* ImageData = new unsigned char[Width * Height * 4];
	for (int k = 0; k < 7; k++)
	{
		int Index = 0;
		char buff[256] = { 0 };
		sprintf_s(buff, 256, "%srendertarget_cylinder%d.jpg", CurrentPath.c_str(), k);
		FIBITMAP* BM = FreeImage_Load(FIF_JPEG, buff);
		if (BM == nullptr)
		{
			printf("Can not load image file...\n");
			exit(1);
		}
		for (int i = 0; i < Height; i++)
		{
			for (int j = 0; j < Width; j++)
			{
				RGBQUAD quad;
				FreeImage_GetPixelColor(BM, j, i, &quad);
				ImageData[Index++] = quad.rgbRed;// ImageBlock::ToRGB(quad.rgbRed, quad.rgbGreen, quad.rgbBlue);
				ImageData[Index++] = quad.rgbGreen;
				ImageData[Index++] = quad.rgbBlue;
				ImageData[Index++] = 255;
			}
		}
		ApplyImageData(k, ImageData);
		FreeImage_Unload(BM);
	}

	SAFE_DELETE_ARRAY(ImageData);

	RenderOneFrame();

	unsigned char* BackBuffer = nullptr;
	FIBITMAP* BM = FreeImage_Allocate(PanoramaWidth, PanoramaHeight, 32);
	BackBuffer = FreeImage_GetBits(BM);

	CopyPanoramaDataFromGPU(BackBuffer);
	char buff[256] = { 0 };
	sprintf_s(buff, 256, "%stest.bmp", CurrentPath.c_str());

	FreeImage_Save(FIF_BMP, BM, buff);
	FreeImage_Unload(BM);

	ReleaseInterface();
}

int main()
{
	printf("Hello, World\n");
	
	StartBuildPanorama();

	return 0;
}