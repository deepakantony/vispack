#include <iostream>
#include <fstream>
#include <utility>
#include <queue>
#include <list>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>

#include <image.h>
#include <imageRGBA.h>
#include <imagefile.h>
#include <rgba.h>

struct ComponentInfo {
	int size;
	int label;
	int labelMaxComponent;
	int labelMaxComponentTemp;
};

const float colorIntensity[] = {1.f, 0.5f, 0.25f, 0.75f, 0.125f, 0.625f, 0.375f, 0.875f};
const int colorIntensitySize = 8;
const rgba colors[] = {rgba(255, 255, 255), rgba(50,100,255), 
	rgba(255, 0, 0), rgba(0, 255, 0), rgba(0,0,255),
	rgba(255, 255, 0), rgba(0, 255, 255), rgba(255, 0, 255)};
const int colorsSize = 8;

using namespace std;

ComponentInfo componentInfo(VISImage<float> labelImg, int label, int labelSize, int sizeThreshold, ComponentInfo* oldInfo) {
	int *labels = (int *)malloc(sizeof(int)*labelSize);
	for(int i = 0; i < labelSize; i++) {
		labels[i] = 0;
	}
	
	int h = labelImg.height();
	int w = labelImg.width();
	int size = 0;
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			if(labelImg(i,j) == label) {
				size++;
				if((i-1) >= 0 && ((int)labelImg.peek(i-1, j)) != label)
					labels[(int)labelImg.peek(i-1, j)]++;
				if((j-1) >= 0 && ((int)labelImg.peek(i, j-1)) != label)
					labels[(int)labelImg.peek(i, j-1)]++;
				if((i+1) < w && ((int)labelImg.peek(i+1, j)) != label)
					labels[(int)labelImg.peek(i+1, j)]++;
				if((j+1) < h && ((int)labelImg.peek(i, j+1)) != label)
					labels[(int)labelImg.peek(i, j+1)]++;
			}
		}
	
	int maxComponent = label;
	for(int i = 0; i < labelSize; i++)
		if(labels[i] > labels[maxComponent] && oldInfo[i].size > sizeThreshold)
			maxComponent = i;
	
	int maxComponentTemp = label;
	if(maxComponent == label) {
		for(int i = 0; i < labelSize; i++)
		if(labels[i] > labels[maxComponentTemp])
			maxComponentTemp = i;
	}
	
	ComponentInfo info;
	info.label = label;
	info.size = size;
	info.labelMaxComponent = maxComponent;
	info.labelMaxComponentTemp = maxComponentTemp;
	return info;
}

int findMaxComponent(ComponentInfo* info, int labelSize, int label, int sizeThreshold) {
	int lbl = label;
	if(label == info[label].labelMaxComponent && info[label].size <= sizeThreshold) {
		while(lbl == info[lbl].labelMaxComponent) {
			int newlbl = info[lbl].labelMaxComponentTemp;
			if(newlbl < 0 || newlbl >= labelSize || newlbl == lbl )
				return label;
			lbl = newlbl;
		}
		return info[lbl].labelMaxComponent;
	}
	else return info[label].labelMaxComponent;
}

ComponentInfo* constructInfo(VISImage<float> labelImg, int labelSize, int sizeThreshold, ComponentInfo* oldInfo) {
	ComponentInfo* infoArray = (ComponentInfo*)malloc(sizeof(ComponentInfo)*labelSize);
	for(int i = 0; i < labelSize; i++) {
		infoArray[i] = componentInfo(labelImg, i, labelSize, sizeThreshold, oldInfo);
	}
	cout << "Initialized the component info" << endl;
	// get the component with the maximum size
	for(int i = 0; i < labelSize; i++) {
		infoArray[i].labelMaxComponent = findMaxComponent(infoArray, labelSize, i, sizeThreshold);
	}
	cout << "Max component for all labels found" << endl;
	return infoArray;
}

VISImage<rgba> createColorImage(VISImage<float> labelImg, float labelSize) {
	VISImage<rgba> colorImg = labelImg.createToSize();
	colorImg = 0;
	
	int h = labelImg.height();
	int w = labelImg.width();
	int ciIdx = 0;
	int cIdx = 0;
	for(float i = 0.f; i < labelSize; i+=1.f) {
		if((((int)i)%8)==0) {
			ciIdx = (ciIdx+1)%colorIntensitySize;
		}
		
		for(int c = 0; c < w; c++)
			for(int r = 0; r < h; r++)
				if(labelImg.peek(c,r) == i) {
					colorImg.poke(c,r) = colors[cIdx]*colorIntensity[ciIdx];
				}
		
		cIdx = (cIdx+1)%colorsSize;
	}
	
	return colorImg;
}

template <class T>
float *histogram(VISImage<T> I, int n, T min, T max) {
	float *H = (float *)malloc(sizeof(float)*n);
	
	for(int i = 0; i < n; i++)
		H[i] = 0.0f;
	
	T mulfact = (n-1)/(max-min);
	int h = I.height();
	int w = I.width();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			T intensity = I.peek(i,j);
			if( intensity >= min && intensity <= max)
				H[(int)((intensity-min)*mulfact)] += 1.0f;
		}
	
	return H;
}

template<class T>
VISImage<byte> threshold(VISImage<T> I, T min, T max) {
	VISImage<byte> outputImg = I.createToSize();
	outputImg = (T)0.0;
	int h = I.height();
	int w = I.width();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			if(I.peek(i,j) >= min && I.peek(i,j) <= max)
				outputImg.poke(i,j) = 255;
	
	return outputImg;
}

template<class T>
int floodfill(VISImage<byte> &I, pair<int,int> seed, int label, VISImage<T> &labelImg) {
	queue< pair<int,int> > Q;
	
	byte cmp = 255;
	
	if(I.peek(seed.first, seed.second) != cmp)
		return 0;
	
	int size = 0;
	Q.push(seed);
	int h = I.height();
	int w = I.width();
	while(!Q.empty()) {
		pair<int,int> front = Q.front();
		if(I.peek(front.first, front.second) == cmp) {
			labelImg.poke(front.first, front.second) = (T)label;
			I.poke(front.first, front.second) = 0;
			size++;
		}
		Q.pop();
		// Check the four neighbours
		if((front.first-1) >= 0 && I.peek(front.first-1, front.second) == cmp) {
			labelImg.poke(front.first-1, front.second) = (T)label;
			I.poke(front.first-1, front.second) = 0;
			size++;
			pair<int,int> p(front.first-1, front.second);
			Q.push(p);
		}
		if((front.first+1) < w && I.peek(front.first+1, front.second) == cmp) {
			labelImg.poke(front.first+1, front.second) = (T)label;
			I.poke(front.first+1, front.second) = 0;
			size++;
			pair<int,int> p(front.first+1, front.second);
			Q.push(p);
		}
		if((front.second-1) >= 0 && I.peek(front.first, front.second-1) == cmp) {
			labelImg.poke(front.first, front.second-1) = (T)label;
			I.poke(front.first, front.second-1) = 0;
			size++;
			pair<int,int> p(front.first, front.second-1);
			Q.push(p);
		}
		if((front.second+1) < h && I.peek(front.first, front.second+1) == cmp) {
			labelImg.poke(front.first, front.second+1) = (T)label;
			I.poke(front.first, front.second+1) = 0;
			size++;
			pair<int,int> p(front.first, front.second+1);
			Q.push(p);
		}
	}
	
	return size;
}

template<class T>
int connectedComponents(VISImage<byte> I, VISImage<T> &labelImg, int initialLabel, list<ComponentInfo> &infoList) {
	int label = initialLabel;
	int h = I.height();
	int w = I.width();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			if(I.peek(i,j) == 255) {
				pair<int, int> seed(i, j);
				int size = floodfill(I, seed, label, labelImg);
				ComponentInfo info;
				info.size = size;
				info.label = label;
				infoList.push_back(info);
				label++;
			}
	
	return label;
}

template<class T>
int connectedComponents(VISImage<byte> I, VISImage<T> &labelImg, int initialLabel) {
	T label = initialLabel;
	int h = I.height();
	int w = I.width();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++)
			if(I.peek(i,j) == 255) {
				pair<int, int> seed(i, j);
				floodfill(I, seed, label, labelImg);
				label++;
			}
	
	return label;
}

int findLabel(ComponentInfo *labelList, int label, int labelListSize) {
	for(int i = 0; i < labelListSize; i++)
		if(labelList[i].label == label) {
			return i;
		}
	return -1;
}

VISImage<float> denoise(VISImage<float> labelImg, int sizeThreshold, ComponentInfo *labelList, int labelSize) {
	VISImage<float> revisedLblImg = labelImg.createToSize();
	
	int h = labelImg.height();
	int w = labelImg.width();
	for(int i = 0; i < w; i++)
		for(int j = 0; j < h; j++) {
			float label = labelImg.peek(i,j);
			int index = findLabel(labelList, label, labelSize);
			//if(labelList[index].size <= sizeThreshold)
				//cout << labelList[index].size << "," << sizeThreshold << endl;
			if(index != -1 && labelList[index].size <= sizeThreshold)
				revisedLblImg.poke(i,j) = labelList[index].labelMaxComponent;
			else
				revisedLblImg.poke(i,j) = label;
		}
	
	return revisedLblImg;
}

void histogramMenu(VISImage<float> I) {
	cout << "Histogram processing..." << endl;
	cout << "Please input your parameters for histogram construction below" << endl;
	cout << "Number of bins: ";
	int n;
	cin >> n;
	cout << "Minimum intensity value: ";
	float min;
	cin >> min;
	cout << "Maximum intensity value: ";
	float max;
	cin >> max;
	float *H = histogram(I, n, min, max);
	cout << "Histogram is constructed... below are the bin values" << endl;
	for(int i = 0; i < (n-1); i++)
		cout << H[i] << ", ";
	cout << H[n-1] << endl;
	
	float sum = 0.0f;
	for(int i = 0; i < n; i++)
		sum += H[i];
	cout << "Sum of the histogram values: " << sum << endl;
	cout << "Number of pixels: " << I.height() * I.width() << endl;
	
	ofstream output;
	output.open("histogram.txt");
	for(int i = 0; i < (n-1); i++)
		output << H[i] << "\n";
	output << H[n-1] << endl;
	output.close();
	
	cout << "Histogram data has been written into the file histogram.txt" << endl;
}

void thresholdMenu(VISImage<float> I) {
	cout << "Dual thresholding..." << endl;
	cout << "Please enter parameters required for dual thresholding" << endl;
	cout << "Minimum intensity (lower limit): ";
	float min;
	cin >> min;
	cout << "Maximum intensity (upper limit): ";
	float max;
	cin >> max;
	VISImage<float> T(threshold(I, min, max));
	VISImageFile imFile;
	imFile.write_jpeg(T, "dual_threshold.jpeg", 100);
	
	cout << "Thresholded file is written into dual_threshold.jpeg" << endl;
}

void connectedCompMenu(VISImage<float> I) {
	cout << "For component analysis you need to threshold the image." << endl;
	cout << "You can enter multiple thresholds, please enter threshold values between min and max." << endl;
	vector<int> thresholdList;
	char ch;
	do {
		int th;
		cout << "Threshold value: ";
		cin >> th;
		thresholdList.push_back(th);
		cout << "Do you want to enter more thresholds (y/n) ?" << endl;
		cin >> ch;
	} while( ch != 'n' && ch != 'N');

	int min = I.min();
	int max = I.max();
	
	thresholdList.push_back(min);
	thresholdList.push_back(max);
	sort(thresholdList.begin(), thresholdList.end());

	string filenameTH("CC_D_TH_");
	string filenameC("CC_final");
	string ext(".jpg");
	VISImageFile imFile;
	int imgN = 0;
	
	// label vars
	VISImage<float> labelImg = I.createToSize();
	labelImg = 0;
	int label = 0;
	
	cout << "Thresholding the image" << endl;
	int sizeTL = thresholdList.size();
	for(int i = 1; i < sizeTL; i++) {
		int th = thresholdList.at(i);
		if(th >= min && th <= max) {
			VISImage<byte> T = threshold(I, (float)min, (float)th);
			min = th;
			
			stringstream ss;
			ss << imgN;
			imgN++;
			string filename = filenameTH + ss.str() + ext;
			imFile.write_jpeg(VISImage<float>(T), filename.c_str(), 100);
			
			label = connectedComponents(T, labelImg, label);
		}
	}
	
	cout << "Analysis complete." << endl;
	
	cout << "Number of components: " << label << endl;
	cout << "writing the labelled image as a colored jpeg..." << endl;
	VISImageRGBA colorImg(createColorImage(labelImg, label));
	imFile.write_jpeg(colorImg, "comp_anlys_labeled.jpg", 100);
	
	cout << "Color image written" << endl;
}

void connectedCompMenu2(VISImage<float> I) {
	cout << "For component analysis you need to threshold the image." << endl;
	cout << "Please enter parameters required for dual thresholding" << endl;
	cout << "Minimum intensity (lower limit): ";
	float min;
	cin >> min;
	cout << "Maximum intensity (upper limit): ";
	float max;
	cin >> max;
	VISImage<byte> T = threshold(I, min, max);
	VISImageFile imFile;
	imFile.write_tiff(T, "comp_anlys_dual_threshold.tif");
	
	cout << "Thresholded file is written into comp_alys_dual_threshold.tif" << endl;
	cout << "Analysing connected components..." << endl;
	
	VISImage<float> labelImg = T.createToSize();
	int label = connectedComponents(T, labelImg, 0);
	VISImage<byte> invT = ((byte)255) - T;
	label = connectedComponents(invT, labelImg, label);
	
	cout << "Analysis complete." << endl;
	
	cout << "Number of components: " << label << endl;
	cout << "writing the labelled image as a colored tiff..." << endl;
	VISImageRGBA colorImg(createColorImage(labelImg, label));
	imFile.write(colorImg, "comp_anlys_labeled.tif");
	
	cout << "Color image written" << endl;
}

void denoiseMenu2(VISImage<float> I) {
	cout << "For topological denoising you need to do the component analysis first" << endl;
	cout << "For component analysis you need to threshold the image." << endl;
	
	string filenameTH("top_denoise_D_TH_");
	string filenameC("top_denoise_CC_");
	string ext(".tif");
	int i = 0;
	char ch = 'n';
	
	cout << "After the 1st thresholded image only foreground component processing is done" << endl;
	cout << "Please select appropriate threshold the first time so that the background is separated properly" << endl;
	bool firstTime = true;
	VISImage<float> labelImg = I.createToSize();
	list<ComponentInfo> infoList;
	int label = 0;
	do {
		cout << "Every new thresholded image will just paste its components into the old image!" << endl;
		cout << "So, after the first thresholding, select thresholds that would create sub-components to the already created components" << endl;
		cout << "Please enter parameters required for dual thresholding" << endl;
		cout << "Minimum intensity (lower limit): ";
		float min;
		cin >> min;
		cout << "Maximum intensity (upper limit): ";
		float max;
		cin >> max;
		VISImage<byte> T = threshold(I, min, max);
		VISImageFile imFile;
		stringstream ss;
		ss << i;
		i++;
		string file = filenameTH + ss.str() + ext;
		imFile.write_tiff(T, file.c_str());
		
		cout << "Thresholded file is written into " << file << endl;
		
		
		cout << "Analysing connected components..." << endl;
	
		
		label = connectedComponents(T, labelImg, label, infoList);
		if(firstTime) {
			VISImage<byte> invT = ((byte)255) - T;
			label = connectedComponents(invT, labelImg, label, infoList);
			firstTime = false;
		}
		
		cout << "Analysis complete." << endl;
		
		cout << "Number of components: " << label << endl;
		cout << "writing the labelled image as a colored tiff..." << endl;
		VISImageRGBA colorImg(createColorImage(labelImg, label));
		file = filenameC + ss.str() + ext;
		imFile.write(colorImg, file.c_str());
		
		cout << "Color image written" << endl;
		cout << "Do you want to add another thresholded image? (y/n)";
		cin >> ch;
	} while(ch != 'n' && ch != 'N');
		
	cout << "Enter the size threshold for the components for denoising: ";
	int size;
	cin >> size;
	
	ComponentInfo *oldInfo = (ComponentInfo*)malloc(sizeof(ComponentInfo)*label);
	i = 0;
	for(list<ComponentInfo>::iterator it = infoList.begin(); it != infoList.end(); it++) {
		oldInfo[i].size = (*it).size;
		oldInfo[i].label = (*it).label;
		i++;
	}

	VISImage<float> revisedLblImg = denoise(labelImg, size, constructInfo(labelImg, label, size, oldInfo), label);
	
	cout << "Writing the revised labelled image as a colored tiff..." << endl;
	VISImageRGBA colorImg(createColorImage(revisedLblImg, label));
	VISImageFile imFile;
	imFile.write(colorImg, "top_denoise_final.tif");
	
	cout << "Final image written" << endl;
}

void denoiseMenu(VISImage<float> I) {
	cout << "For topological denoising you need to do the component analysis first" << endl;
	
	cout << "For component analysis you need to threshold the image." << endl;
	cout << "You can enter multiple thresholds, please enter threshold values between min and max." << endl;
	vector<int> thresholdList;
	char ch;
	do {
		int th;
		cout << "Threshold value: ";
		cin >> th;
		thresholdList.push_back(th);
		cout << "Do you want to enter more thresholds (y/n) ?" << endl;
		cin >> ch;
	} while( ch != 'n' && ch != 'N');

	int min = I.min();
	int max = I.max();
	
	thresholdList.push_back(min);
	thresholdList.push_back(max);
	sort(thresholdList.begin(), thresholdList.end());

	string filenameTH("TD_D_TH_");
	string filenameC("TD_final");
	string ext(".jpg");
	VISImageFile imFile;
	int imgN = 0;
	
	// label vars
	VISImage<float> labelImg = I.createToSize();
	labelImg = 0;
	int label = 0;
	list<ComponentInfo> infoList;
	
	cout << "Thresholding the image" << endl;
	int sizeTL = thresholdList.size();
	for(int i = 1; i < sizeTL; i++) {
		int th = thresholdList.at(i);
		if(th >= min && th <= max) {
			VISImage<byte> T = threshold(I, (float)min, (float)th);
			min = th;
			
			stringstream ss;
			ss << imgN;
			imgN++;
			string filename = filenameTH + ss.str() + ext;
			imFile.write_jpeg(VISImage<float>(T), filename.c_str(), 100);
			
			label = connectedComponents(T, labelImg, label, infoList);
		}
	}
	
	cout << "Connected component analysis complete." << endl;
	
	cout << "Number of components: " << label << endl;
	cout << "writing the labelled image as a colored jpg..." << endl;
	VISImageRGBA colorImg(createColorImage(labelImg, label));
	imFile.write_jpeg(colorImg, "TD_label_cc.jpg", 100);
	
	cout << "Components image written" << endl;
		
	cout << "Enter the size threshold for the components for denoising: ";
	int size;
	cin >> size;
	
	cout << "Denoising the components" << endl;
	ComponentInfo *oldInfo = (ComponentInfo*)malloc(sizeof(ComponentInfo)*label);
	int i = 0;
	for(list<ComponentInfo>::iterator it = infoList.begin(); it != infoList.end(); it++) {
		oldInfo[i].size = (*it).size;
		oldInfo[i].label = (*it).label;
		i++;
	}
	
	ComponentInfo *newInfo = constructInfo(labelImg, label, size, oldInfo);
	
	cout << "Component Info calculated" << endl;

	VISImage<float> revisedLblImg = denoise(labelImg, size, newInfo, label);
	cout << "Found the revised label image" << endl;
	
	cout << "Writing the revised labelled image as a colored jpeg..." << endl;
	VISImageRGBA colorImg2(createColorImage(revisedLblImg, label));
	imFile.write_jpeg(colorImg2, "TD_final.jpg", 100);
	
	cout << "Final image written" << endl;
}

int main(int argc, char** argv) {
	VISIm im;
	VISImageFile imgFile;
	VISImage<float> I;
	if(argc > 1 && (im = imgFile.read(argv[1])).isValid()) {
		I = VISImage<float>(im);
	}
	else {
		cout << "Image file could not be read or no file was specified!" << endl;
		return -1;
	}
	cout << "Image's size: " << I.height() << "x" << I.width() << endl;
	cout << "Image's max intensity is " << I.max() << " and min intensity is " << I.min() << endl;
	
	while(true) {
		cout << "Here are your options... Please enter the corresponding number." << endl;
		
		cout << "1. Histogram" << endl;
		cout << "2. Dual Thresholding" << endl;
		cout << "3. Connected components" << endl;
		cout << "4. Topological denoising" << endl;
		cout << "5. Exit" << endl;
		
		int option;
		cin >> option;
		
		switch(option) {
			case 1: histogramMenu(I); break;
			case 2: thresholdMenu(I); break;
			case 3: connectedCompMenu(I); break;
			case 4: denoiseMenu(I); break;
			case 5: return 0;
			default: 
				cout << "Wrong option entered!" << endl;
				break;
		}
	}
	return 0;
}
