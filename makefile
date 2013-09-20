all : img.exe convolve.dll

img.exe : img.cs ImageUtils.cs
	csc /unsafe /checked- /O img.cs ImageUtils.cs

convolve.dll : convolve.obj
	link /dll convolve.obj

convolve.obj : convolve.cpp
	cl /c /Ob2 /O2 /Oi /fp:fast convolve.cpp

clean :
	-del img.exe
	-del convolve.dll
	-del convolve.obj
	-del convolve.exp
	-del convolve.lib
