/***********************************************
 * File: Interface.h
 *
 * Author: LYG
 * Date: ÆßÔÂ 2018
 *
 * Purpose:
 *
 * 
 **********************************************/
#pragma once
#define				API_EXPORT				extern "C" __declspec(dllexport)
#define				DLL_CALLCONV			__stdcall

enum Renderer
{
	RendererCPU,
	RendererD3D11,
	RendererD3D12,
};
// Set the renderer type
API_EXPORT void DLL_CALLCONV SetPanoramaRendererType(Renderer Type);
// Set the panorama image width and height you want to build
API_EXPORT void DLL_CALLCONV SetPanoramaWidthAndHeight(int Width, int Height);
// Set the input image width and height
API_EXPORT void DLL_CALLCONV SetInputImageWidthAndHeight(int Width, int Height);
// Set the input image projection parameters, all the input image has the same.
API_EXPORT void DLL_CALLCONV SetInputImageProjectionParameters(float Fovy, float NearClip, float FarClip);
// Set the count of the input images
API_EXPORT void DLL_CALLCONV SetInputImageCount(int ImageCount);
// Set the input image rotation by index.it's a quaternion, in x,y,z,w order.
API_EXPORT void DLL_CALLCONV SetInputImageRotation(int Index, float* InputImageQuaternionXYZW);
// Initialise the context and environment such as d3d11/d3d12 etc...
API_EXPORT bool DLL_CALLCONV Initialise();
// Set the input image data by index
// ImageData must be r8g8b8a8 format, the buffer length must be width * height * 4.
API_EXPORT void DLL_CALLCONV ApplyImageData(int Index, unsigned char* ImageData);
// Render one frame.
API_EXPORT void DLL_CALLCONV RenderOneFrame();
// Get the output panorama image data. the buffer length must be 4 * PanoramaWidth * PanoramaHeight
API_EXPORT void DLL_CALLCONV CopyPanoramaDataFromGPU(unsigned char* Buffer);
API_EXPORT void DLL_CALLCONV ReleaseInterface();
/*
usage:
int main()
{
	SetPanoramaRendererType();
	SetPanoramaWidthAndHeight();
	SetInputImageWidthAndHeight();
	SetInputImageProjectionParameters();
	SetInputImageCount();
	SetInputImageRotation();
	Initialise();

	ApplyImageData();
	RenderOneFrame();
	CopyPanoramaDataFromGPU();
	// repeat
	ApplyImageData();
	RenderOneFrame();
	CopyPanoramaDataFromGPU();
	.......

	ReleaseInterface();
	return 0;
}	

*/