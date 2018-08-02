PanoramaPlugin
=====
There is a panorama plugin, which can build panorama image using all kind of resolution of input images.
------
~~~cpp
#pragma once
#define				API_EXPORT			extern "C" __declspec(dllexport)
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

// Set the input image projection parameters, all the input images are the same.
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
~~~
In the most render engines, building a panorama image in the following step:
1. build the camera cubemap image, which in six faces, each face is a quad, may be in the screen resolution of 1024 * 1024.
2. sample the cube map image to panorama image.
but in this way, the cubemap render must be in a render target, not screen capture. and in the most render engine like unreal 4, may be lost many post process volume effects.

in my plugin, i can build panorama image in all kind of screen resolution, so you can capture screen to build the panorama, without lose any post process volume effects.

any one who like this plugin, please contact me, not for free. in this folder, you can use the debug dll to do any test.
my email address: lygyue@126.com

~Any one who want to compile the application, please config the platform to debug + x86 mode.~
~The best renderer is d3d11. d3d12 need win10, and cpu using 4 threads to do the sample, it's slower than d3d11.~
~Some one may think the plugin run slowly, in fact ,it's very fast.the initialise may be a bit slow, it run once at the begining.~
~I don't do enough tests, any issues please contact me, thankyou!~
