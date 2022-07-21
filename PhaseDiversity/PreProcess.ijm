run("Close All");
pathIn = "D:/Data/20220718_Test/PD_Batch_0p03_step0p5um/"
pathOut = pathIn + "processed/"
File.makeDirectory(pathOut);
bgName = "AVG_bg.tif";
fileBg = "D:/Data/20220718_Test/bg_g0_30ms/" + bgName;
print("Processing start ...");
for (i = 1; i <= 52; i++) {
	print("Processing img #: " + i);
	fileName = "A_1.tif";
	fileIn = pathIn + "A_" + i +"/" + fileName;
	fileOut = pathOut + "A_" + i +"_1.tif";
	open(fileIn);
	open(fileBg);
	imageCalculator("Subtract create stack", fileName, bgName);
	makeRectangle(374,435, 128, 128);
	run("Crop");
	saveAs("Tiff", fileOut);
	run("Close All");
	
	fileName = "A_-1.tif";
	fileIn = pathIn + "A_" + i +"/" + fileName;
	fileOut = pathOut + "A_" + i +"_-1.tif";
	open(fileIn);
	open(fileBg);
	imageCalculator("Subtract create stack", fileName, bgName);
	makeRectangle(374, 435, 128, 128);
	run("Crop");
	saveAs("Tiff", fileOut);
	run("Close All");
}
print("Processing completed !!!");

