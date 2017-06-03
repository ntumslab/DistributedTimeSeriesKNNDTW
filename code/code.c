#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <float.h>

//#define EQUAL_SIZE
//#define ORACLE_INITIALIZATION
//#define RANDOM_INITIALIZATION
//#define COMPRESSION

typedef struct _FileStream{
	double *stream;
	int *elementSizes;
	int size;
} FileStream;

FileStream* fileStreamConstruct(){
	FileStream *fileStream = (FileStream*)malloc(sizeof(FileStream));
	fileStream->stream = (double*)malloc(sizeof(double));
	fileStream->elementSizes = (int*)malloc(sizeof(int));
	fileStream->size = 0;
	return fileStream;
}

void fileStreamDestruct(FileStream *fileStream){
	if(fileStream == NULL){
		return;
	}
	free(fileStream->stream);
	free(fileStream->elementSizes);
	free(fileStream);
}

int countFileSize(char *fileName){
	FILE *file = fopen(fileName, "r");
	fseek(file, 0, SEEK_END);
	int size = ftell(file);
	fclose(file);
	return size;
}

void addValueToStream(unsigned char *stream, int *pointer, int value, int valueSize){
	for(int d = valueSize - 1; d >= 0; d --, *pointer += 1){
		int ptrByte = *pointer / 8;
		int ptrBit = *pointer % 8;

		unsigned char bit = (unsigned char)((value >> d) & 1);
		unsigned char newBit = bit << (7 - ptrBit);
		stream[ptrByte] |= newBit;
	}
}

char *ORIGINAL_FILE_NAME;
char *COMPRESSED_FILE_NAME;

double fileStreamTestCompressionRate(FileStream *fileStream, int *originalSize, int *compressedSize){
	char *originalFileName = ORIGINAL_FILE_NAME;
	char *compressedFileName = COMPRESSED_FILE_NAME;
	
	int totalSize = 0;
	for(int i = 0; i < fileStream->size; i ++){
		totalSize += fileStream->elementSizes[i];
	}
	int bitSize = (totalSize + 7) / 8;
	int pointer = 0;
	unsigned char *stream = (unsigned char*)malloc(sizeof(unsigned char) * bitSize);
	for(int b = 0; b < bitSize; b ++){
		stream[b] = 0;
	}

	for(int i = 0; i < fileStream->size; i ++){
		addValueToStream(stream, &pointer, (int)fileStream->stream[i], fileStream->elementSizes[i]);
	}

	FILE *file = fopen(originalFileName, "wb");
	fwrite(stream, sizeof(unsigned char), bitSize, file);
	fclose(file);
	free(stream);

	*originalSize = countFileSize(originalFileName);

	char compressionCommand[200];
	sprintf(compressionCommand, "zip -q -r %s %s", compressedFileName, originalFileName);
	system(compressionCommand);

	*compressedSize = countFileSize(compressedFileName);

	return (double)*compressedSize / *originalSize;
}

//================

typedef enum{
	LB_MS,
	LB_KEOGH,
	LB_KIM
} DTWLowerBound;

typedef struct{
	double min;
	double max;
	int size;
} Segment;

typedef struct _Sensor{
	int index;
	double label;
	int realSegmentCount;
	int segmentCount;
	Segment *realSegments;
	Segment *segments;
	bool isCandidate;
	struct _Sensor *coarseSensor;
} Sensor;

typedef struct{
	int queryLength;
	int candidateLength;
	Sensor *query;
	Sensor *candidate;
	int *lowerConstraints;
	int *upperConstraints;
	double constraintHalfWidth;
} DTW;

typedef struct{
	int queryLength;
	Sensor *query;
	Sensor *candidate;
	double *lowerConstraints;
	double *upperConstraints;
	double constraintHalfWidth;
} LBKeogh;

typedef struct{
	Sensor *query;
	Sensor *candidate;
} LBKim;

typedef struct{
	int index;
	int sensorCount;
	Sensor **sensors;
	int candidateCount;
} Site;

typedef struct{
	int realSensorCount;
	int sensorCount;
	int siteCount;
	int queryLength;
	Sensor **sensors;
	Site **sites;
	Sensor *query;
} Dataset;

typedef struct{
	int sensorIndex;
	int siteIndex;
	double dtw;
	double lowerBound;
	double upperBound;
} KNNData;

typedef struct _Model{
	int K;
	int bandwidth;
	int dataUnitBits;
	int firstLevelSegmentCount;
	int subSegmentCount;
	KNNData *knnList;
	Dataset *dataset;
	void (*knnStart)(struct _Model*, int);
	FileStream *fileStream;
} Model;

//================

DTW* DTWConstruct(Sensor*, Sensor*);
void DTWDestruct(DTW*);

Sensor* sensorConstruct(int, int);
void sensorDestruct(Sensor*);

Site* siteConstruct();
void siteDestruct(Site*);

Dataset* datasetConstruct(char*);
void datasetDestruct(Dataset*);

void modelBandwidthIncrement(Model*, int, bool);

//================

// Computes distance between two data points.
double boundDistance(Segment *a, Segment *b){
	double distance;
	
	if(a->min > b->max){
		distance = (a->min - b->max) * (a->min - b->max);
	}
	else if(b->min > a->max){
		distance = (b->min - a->max) * (b->min - a->max);
	}
	else{
		distance = 0.0;
	}
	distance *= (a->size < b->size)? a->size: b->size;
	
	return distance;
}

//================
// LB_KimFL

// Runs LB_KimFL.
double LBKimStart(LBKim *lbkim){
	double firstDistance = boundDistance(&(lbkim->query->segments[0]), &(lbkim->candidate->segments[0]));
	double lastDistance = boundDistance(&(lbkim->query->segments[lbkim->query->segmentCount - 1]), &(lbkim->candidate->segments[lbkim->candidate->segmentCount - 1]));
	return (firstDistance > lastDistance) ? firstDistance : lastDistance;
}

LBKim* LBKimConstruct(Sensor *query, Sensor *candidate){
	LBKim *lbkim = (LBKim*)malloc(sizeof(LBKim));
	lbkim->query = query;
	lbkim->candidate = candidate;
	return lbkim;
}

void LBKimDestruct(LBKim *lbkim){
	free(lbkim);
}

//================
// LB_Keogh

// Runs LB_Keogh.
// threshold: For early abandoning.
double LBKeoghStart(LBKeogh *lbkeogh, double threshold){
	double lowerBound = 0.0;
	Segment queryBounds;
	queryBounds.size = 1;
	for(int t = 0; t < lbkeogh->queryLength; t ++){
		queryBounds.min = lbkeogh->lowerConstraints[t];
		queryBounds.max = lbkeogh->upperConstraints[t];
		lowerBound += boundDistance(&queryBounds, &lbkeogh->candidate->segments[t]);
		if(lowerBound > threshold){
			return lowerBound; // Early abandoning: Returns any value greater than the threshold.
		}
	}
	return lowerBound;
}

LBKeogh* LBKeoghConstruct(Sensor *query, Sensor *candidate){
	LBKeogh *lbkeogh = (LBKeogh*)malloc(sizeof(LBKeogh));
	lbkeogh->query = query;
	lbkeogh->candidate = candidate;
	lbkeogh->queryLength = query->segmentCount;
	lbkeogh->constraintHalfWidth = 0.025;
	
	lbkeogh->lowerConstraints = (double*)malloc(sizeof(double) * lbkeogh->queryLength);
	lbkeogh->upperConstraints = (double*)malloc(sizeof(double) * lbkeogh->queryLength);
	for(int t = 0; t < lbkeogh->queryLength; t ++){
		int leftBorder = t - (int)(lbkeogh->queryLength * lbkeogh->constraintHalfWidth);
		if(leftBorder < 0){
			leftBorder = 0;
		}
		
		int rightBorder = t + (int)(lbkeogh->queryLength * lbkeogh->constraintHalfWidth) + 1;
		if(rightBorder > lbkeogh->queryLength){
			rightBorder = lbkeogh->queryLength;
		}
		
		double min, max;
		bool minNull = true, maxNull = true;
		for(int q = leftBorder; q < rightBorder; q ++){
			if(minNull || min > lbkeogh->query->segments[q].min){
				min = lbkeogh->query->segments[q].min;
				minNull = false;
			}
			if(maxNull || max < lbkeogh->query->segments[q].max){
				max = lbkeogh->query->segments[q].max;
				maxNull = false;
			}
		}
		
		lbkeogh->lowerConstraints[t] = min;
		lbkeogh->upperConstraints[t] = max;
	}
	
	return lbkeogh;
}

void LBKeoghDestruct(LBKeogh *lbkeogh){
	free(lbkeogh->lowerConstraints);
	free(lbkeogh->upperConstraints);
	free(lbkeogh);
	lbkeogh = NULL;
}

//================
// LB_MS or normal DTW
// When every segment has max == min and len == 1, LB_MS outputs are the same as DTW.

bool DTWIsLegalCell(DTW *dtw, int q, int c){
	return q >= 0 && q < dtw->queryLength && c >= dtw->lowerConstraints[q] && c < dtw->upperConstraints[q];
}

// Runs LB_MS.
// threshold: For early abandoning.
double DTWStart(DTW *dtw, double threshold){
	double warpingTable[2][dtw->candidateLength];
	for(int q = 0; q < dtw->queryLength; q ++){
		bool exceedThreshold = true;
		for(int c = dtw->lowerConstraints[q]; c < dtw->upperConstraints[q]; c ++){
				bool isNull = true;
				double minValue = 0.0;
				if(DTWIsLegalCell(dtw, q - 1, c - 1)){
					minValue = warpingTable[(q - 1) & 1][c - 1];
					isNull = false;
				}
				if(DTWIsLegalCell(dtw, q - 1, c) && (isNull || minValue > warpingTable[(q - 1) & 1][c])){
					minValue = warpingTable[(q - 1) & 1][c];
					isNull = false;
				}
				if(DTWIsLegalCell(dtw, q, c - 1) && (isNull || minValue > warpingTable[q & 1][c - 1])){
					minValue = warpingTable[q & 1][c - 1];
					isNull = false;
				}
				warpingTable[q & 1][c] = boundDistance(&dtw->query->segments[q], &dtw->candidate->segments[c]) + minValue;

				if(warpingTable[q & 1][c] <= threshold){
					exceedThreshold = false;
				}
		}
		if(exceedThreshold){
			return threshold + 1; // Early abandoning: Returns any value greater than the threshold.
		}
	}
	return warpingTable[(dtw->queryLength - 1) & 1][dtw->candidateLength - 1];
}

DTW* DTWConstruct(Sensor *query, Sensor *candidate){
	DTW *dtw = (DTW*)malloc(sizeof(DTW));
	dtw->query = query;
	dtw->candidate = candidate;
	dtw->queryLength = query->segmentCount;
	dtw->candidateLength = candidate->segmentCount;
	dtw->constraintHalfWidth = 0.025; // Sakoe-Chiba constraints half width = 2.5%.
	
	dtw->lowerConstraints = (int*)malloc(sizeof(int) * dtw->queryLength);
	dtw->upperConstraints = (int*)malloc(sizeof(int) * dtw->queryLength);
	double lengthRatio = (double)dtw->candidateLength / dtw->queryLength;
	for(int q = 0; q < query->segmentCount; q ++){
		int leftBound = q - (int)(lengthRatio * dtw->constraintHalfWidth);
		if(leftBound < 0){
			leftBound = 0;
		}

		int rightBound = q + (int)(lengthRatio * dtw->constraintHalfWidth) + 1;
		if(rightBound > dtw->candidateLength){
			rightBound = dtw->candidateLength;
		}

		dtw->lowerConstraints[q] = leftBound;
		dtw->upperConstraints[q] = rightBound;
	}
	
	return dtw;
}

void DTWDestruct(DTW *dtw){
	free(dtw->lowerConstraints);
	free(dtw->upperConstraints);
	free(dtw);
	dtw = NULL;
}

//================
// Sensor-related functions

// Finds the maximum or the minimum in a segment.
// isMin: true => minimum; false => maximum.
double segmentExtremeValue(Segment* segments, int n, bool isMin){
	double extremeValue;
	bool isNull = true;
	if(isMin){
		for(int t = 0; t < n; t ++){
			if(isNull || extremeValue > segments[t].min){
				extremeValue = segments[t].min;
				isNull = false;
			}
		}
	}
	else{
		for(int t = 0; t < n; t ++){
			if(isNull || extremeValue < segments[t].max){
				extremeValue = segments[t].max;
				isNull = false;
			}
		}
	}
	return extremeValue;
}

// Z-normalization: x' = (x - mean) / stddev.
void sensorNormalizeSegments(Sensor *sensor){
	double minMean = 0, maxMean = 0;
	int length = 0;
	for(int t = 0; t < sensor->segmentCount; t ++){
		minMean += sensor->segments[t].min * sensor->segments[t].size;
		maxMean += sensor->segments[t].max * sensor->segments[t].size;
		length += sensor->segments[t].size;
	}
	minMean /= length;
	maxMean /= length;

	double minStd = 0, maxStd = 0;
	for(int t = 0; t < sensor->segmentCount; t ++){
		double minDiff = sensor->segments[t].min - minMean;
		double maxDiff = sensor->segments[t].max - maxMean;
		minStd += minDiff * minDiff * sensor->segments[t].size;
		maxStd += maxDiff * maxDiff * sensor->segments[t].size;
	}
	minStd = sqrt(minStd / length);
	maxStd = sqrt(maxStd / length);

	for(int t = 0; t < sensor->segmentCount; t ++){
		sensor->segments[t].min = (sensor->segments[t].min - minMean) / minStd;
		sensor->segments[t].max = (sensor->segments[t].max - maxMean) / maxStd;
	}
}

// Returns the level-1 coarse query.
// segmentDistribution: Assigns the length (size) of each segment.
// isEqual: true => equal-size segmentation; false => unequal-size segmentation.
Sensor* sensorSegment(Sensor *sensor, int *segmentDistribution, bool isEqual){
	Sensor *newSensor;
	if(isEqual){
		int segmentSize = *segmentDistribution;
		int segmentCount = (sensor->segmentCount + segmentSize - 1) / segmentSize;
		newSensor = sensorConstruct(sensor->index, segmentCount);
		
		for(int t = 0; t < sensor->segmentCount; t += segmentSize){
			Segment *subSegment = sensor->segments + t;
			int thisSegmentSize = (sensor->segmentCount - t >= segmentSize) ? segmentSize : sensor->segmentCount - t;

			double min = segmentExtremeValue(subSegment, thisSegmentSize, true);
			double max = segmentExtremeValue(subSegment, thisSegmentSize, false);
			newSensor->segments[t / segmentSize].min = min;
			newSensor->segments[t / segmentSize].max = max;
			newSensor->segments[t / segmentSize].size = thisSegmentSize;
		}
	}
	else{
		// The terminal value 0 in segmentDistribution means that the remaining data points belong to the last segment.
		int segmentCount = 0;
		for(int segmentSize = *segmentDistribution, t = 0; (segmentSize = segmentDistribution[t]) != 0; t ++){
			segmentCount ++;
		}
		segmentCount ++;
		newSensor = sensorConstruct(sensor->index, segmentCount);

		int u, t, segmentSize;
		for(segmentSize = *segmentDistribution, t = 0, u = 0; (segmentSize = segmentDistribution[t]) != 0; t ++, u += segmentSize){
			Segment *subSegment = sensor->segments + u;
			double min = segmentExtremeValue(subSegment, segmentSize, true);
			double max = segmentExtremeValue(subSegment, segmentSize, false);
			newSensor->segments[t].min = min;
			newSensor->segments[t].max = max;
			newSensor->segments[t].size = segmentSize;
		}
		// Last segment
		segmentSize = sensor->segmentCount - u;
		Segment *subSegment = sensor->segments + u;
		double min = segmentExtremeValue(subSegment, segmentSize, true);
		double max = segmentExtremeValue(subSegment, segmentSize, false);
		newSensor->segments[t].min = min;
		newSensor->segments[t].max = max;
		newSensor->segments[t].size = segmentSize;
	}
	return newSensor;
}

// Analyzes how to cut a time series into unequal-size segments with a simple bottom-up algorithm.
// sensor: A finest time series (max == min, len == 1).
// wantSegmentCount: The parameter N. How many segments the time series is cut into.
int* sensorAnalyzeSegmentBottomUp(Sensor *sensor, int wantSegmentCount){
	int *segmentList = (int*)malloc(sizeof(int));
	int listCount = 0;

	// Initialization
	Segment *analyzedSegments = (Segment*)malloc(sizeof(Segment) * sensor->segmentCount);
	int analyzedSegmentCount = sensor->segmentCount;
	for(int t = 0; t < sensor->segmentCount; t ++){
		analyzedSegments[t] = sensor->segments[t];
	}
	
	while(analyzedSegmentCount > wantSegmentCount){
		double upperBound, lowerBound;
		double minMergeError;
		int argminMergeError;
		for(int t = 0, q = 0; t < analyzedSegmentCount - 1; q += analyzedSegments[t].size, t ++){
			double maxBound = (analyzedSegments[t].max > analyzedSegments[t + 1].max) ? analyzedSegments[t].max : analyzedSegments[t + 1].max;
			double minBound = (analyzedSegments[t].min < analyzedSegments[t + 1].min) ? analyzedSegments[t].min : analyzedSegments[t + 1].min;
			double thisMergeError = 0;
			
			for(int k = q; k < q + analyzedSegments[t].size + analyzedSegments[t + 1].size; k ++){
				double maxDiff = sensor->segments[k].max - maxBound;
				double minDiff = sensor->segments[k].min - minBound;
				thisMergeError += maxDiff * maxDiff + minDiff * minDiff;
			}
			thisMergeError /= analyzedSegments[t].size + analyzedSegments[t + 1].size;
			
			if(t == 0 || minMergeError > thisMergeError){
				minMergeError = thisMergeError;
				argminMergeError = t;
				upperBound = maxBound;
				lowerBound = minBound;
			}
		}
		
		// Merges of two adjacent segments.
		analyzedSegments[argminMergeError].size += analyzedSegments[argminMergeError + 1].size;
		analyzedSegments[argminMergeError].max = upperBound;
		analyzedSegments[argminMergeError].min = lowerBound;
		
		for(int t = argminMergeError + 1; t < analyzedSegmentCount - 1; t ++){
			analyzedSegments[t] = analyzedSegments[t + 1];
		}
		analyzedSegmentCount --;
	}
	
	for(int t = 0; t < analyzedSegmentCount - 1; t ++){
		segmentList = (int*)realloc(segmentList, sizeof(int) * (listCount + 1));
		segmentList[listCount] = analyzedSegments[t].size;
		listCount ++;
	}
	segmentList = (int*)realloc(segmentList, sizeof(int) * (listCount + 1));
	segmentList[listCount] = 0;
	listCount ++;

	free(analyzedSegments);
	return segmentList;
}

// Encodes a data stream.
// sensor: A coarse time series at certain resolution level.
// accurateSensor: Finest query (max == min, len == 1) as the encoding reference.
// subSegmentCount: The parameter R. How many segments at the next level a segment at this level is cut into.
// dataUnitBits: The parameter B. The space of a data point.
double* sensorEncodeSplitInfo(Sensor *sensor, Sensor *accurateSensor, int subSegmentCount, int dataUnitBits, int *bandwidth, FileStream* encodedStream){
	double *transmittedData = (double*)malloc(sizeof(double));
	int dataCount = 0;
	int signalBits = (int)ceil(log(subSegmentCount) / log(2));
	*bandwidth = 0;
	
	for(int t = 0, u = 0; t < sensor->segmentCount; u += sensor->segments[t].size, t ++){
		// General rule 1: Does nothing.
		// General rule 2
		if(sensor->segments[t].size == 2){
			// 0 means q_u < q_{u+1}; 1 otherwise.
			double signal = (accurateSensor->segments[u].min == sensor->segments[t].min) ? 0.0 : 1.0;
			transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 1));
			transmittedData[dataCount] = signal;
			dataCount ++;
			////
			*bandwidth += 1; // 1-bit signal to indicate the order of the min and the max
			////
#ifdef COMPRESSION
			encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 1));
			encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 1));
			encodedStream->stream[encodedStream->size] = signal;
			encodedStream->elementSizes[encodedStream->size] = 1;
			encodedStream->size += 1;
#endif
		}
		// General rule 3
		else if(sensor->segments[t].size > 2 && sensor->segments[t].size <= subSegmentCount){
			int oldMinPosition, oldMaxPosition;
			for(int i = 0; i < sensor->segments[t].size; i ++){
				if(accurateSensor->segments[u + i].min == sensor->segments[t].min){
					oldMinPosition = i;
				}
				if(accurateSensor->segments[u + i].max == sensor->segments[t].max){
					oldMaxPosition = i;
				}
			}
			transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 2));
			transmittedData[dataCount] = (double)oldMinPosition;
			transmittedData[dataCount + 1] = (double)oldMaxPosition;
			dataCount += 2;
			////
			int smallSegmentSignalBits = (int)ceil(log(sensor->segments[t].size) / log(2));
			*bandwidth += 2 * smallSegmentSignalBits; // 2 signals to indicate new locations of max, min.
			////
#ifdef COMPRESSION
			encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 2));
			encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 2));
			encodedStream->stream[encodedStream->size] = oldMinPosition;
			encodedStream->stream[encodedStream->size + 1] = oldMaxPosition;
			encodedStream->elementSizes[encodedStream->size] = smallSegmentSignalBits;
			encodedStream->elementSizes[encodedStream->size + 1] = smallSegmentSignalBits;
			encodedStream->size += 2;
#endif
			for(int i = 0; i < sensor->segments[t].size; i ++){
				if(i == oldMinPosition || i == oldMaxPosition){
					continue;
				}
				transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 1));
				transmittedData[dataCount] = accurateSensor->segments[u + i].min;
				dataCount ++;
				////
				*bandwidth += 1 * dataUnitBits; // extra min / max data
				////
#ifdef COMPRESSION
				encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 1));
				encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 1));
				encodedStream->stream[encodedStream->size] = accurateSensor->segments[u + i].min;
				encodedStream->elementSizes[encodedStream->size] = dataUnitBits;
				encodedStream->size += 1;
#endif
			}
		}
		// General rule 4
		else if(sensor->segments[t].size > subSegmentCount && sensor->segments[t].size < subSegmentCount * 2){
			int oldMinPosition, oldMaxPosition;
			int remainingSize = sensor->segments[t].size - (subSegmentCount - 1);
			double endSubMin = segmentExtremeValue(accurateSensor->segments + (u + (subSegmentCount - 1)), remainingSize, true);
			double endSubMax = segmentExtremeValue(accurateSensor->segments + (u + (subSegmentCount - 1)), remainingSize, false);
			for(int i = 0; i < subSegmentCount - 1; i ++){
				if(accurateSensor->segments[u + i].min == sensor->segments[t].min){
					oldMinPosition = i;
				}
				if(accurateSensor->segments[u + i].max == sensor->segments[t].max){
					oldMaxPosition = i;
				}
			}
			if(endSubMin == sensor->segments[t].min){
				oldMinPosition = subSegmentCount - 1;
			}
			if(endSubMax == sensor->segments[t].max){
				oldMaxPosition = subSegmentCount - 1;
			}
			transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 2));
			transmittedData[dataCount] = (double)oldMinPosition;
			transmittedData[dataCount + 1] = (double)oldMaxPosition;
			dataCount += 2;
			////
			*bandwidth += 2 * signalBits; // 2 signals to indicate new locations of max, min.
			////
#ifdef COMPRESSION
			encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 2));
			encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 2));
			encodedStream->stream[encodedStream->size] = oldMinPosition;
			encodedStream->stream[encodedStream->size + 1] = oldMaxPosition;
			encodedStream->elementSizes[encodedStream->size] = signalBits;
			encodedStream->elementSizes[encodedStream->size + 1] = signalBits;
			encodedStream->size += 2;
#endif
			for(int i = 0; i < subSegmentCount - 1; i ++){
				if(i == oldMinPosition || i == oldMaxPosition){
					continue;
				}
				transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 1));
				transmittedData[dataCount] = accurateSensor->segments[u + i].min;
				dataCount ++;
				////
				*bandwidth += 1 * dataUnitBits; // Extra max, min data points.
				////
#ifdef COMPRESSION
				encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 1));
				encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 1));
				encodedStream->stream[encodedStream->size] = accurateSensor->segments[u + i].min;
				encodedStream->elementSizes[encodedStream->size] = dataUnitBits;
				encodedStream->size += 1;
#endif
			}
			if(subSegmentCount - 1 != oldMinPosition){
				transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 1));
				transmittedData[dataCount] = endSubMin;
				dataCount ++;
				////
				*bandwidth += 1 * dataUnitBits; // Extra max, min data points.
				////
#ifdef COMPRESSION
				encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 1));
				encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 1));
				encodedStream->stream[encodedStream->size] = endSubMin;
				encodedStream->elementSizes[encodedStream->size] = dataUnitBits;
				encodedStream->size += 1;
#endif
			}
			if(subSegmentCount - 1 != oldMaxPosition){
				transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 1));
				transmittedData[dataCount] = endSubMax;
				dataCount ++;
				////
				*bandwidth += 1 * dataUnitBits; // Extra max, min data points.
				////
#ifdef COMPRESSION
				encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 1));
				encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 1));
				encodedStream->stream[encodedStream->size] = endSubMax;
				encodedStream->elementSizes[encodedStream->size] = dataUnitBits;
				encodedStream->size += 1;
#endif
			}	
		}
		// General rule 5
		else if(sensor->segments[t].size >= subSegmentCount * 2){
			double subMin[subSegmentCount];
			double subMax[subSegmentCount];
			int subSize = sensor->segments[t].size / subSegmentCount;
			int remainingSize = sensor->segments[t].size - subSize * (subSegmentCount - 1);
			int oldMinPosition, oldMaxPosition;
			for(int i = 0; i < subSegmentCount; i ++){
				if(i < subSegmentCount - 1){
					subMin[i] = segmentExtremeValue(accurateSensor->segments + (u + subSize * i), subSize, true);
					subMax[i] = segmentExtremeValue(accurateSensor->segments + (u + subSize * i), subSize, false);
				}
				else{
					subMin[i] = segmentExtremeValue(accurateSensor->segments + (u + subSize * i), remainingSize, true);
					subMax[i] = segmentExtremeValue(accurateSensor->segments + (u + subSize * i), remainingSize, false);
				}
				if(subMin[i] == sensor->segments[t].min){
					oldMinPosition = i;
				}
				if(subMax[i] == sensor->segments[t].max){
					oldMaxPosition = i;
				}
			}
			transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 2));
			transmittedData[dataCount] = (double)oldMinPosition;
			transmittedData[dataCount + 1] = (double)oldMaxPosition;
			dataCount += 2;
			////
			*bandwidth += 2 * signalBits; // 2 signals to indicate new locations of max, min.
			////
#ifdef COMPRESSION
			encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 2));
			encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 2));
			encodedStream->stream[encodedStream->size] = oldMinPosition;
			encodedStream->stream[encodedStream->size + 1] = oldMaxPosition;
			encodedStream->elementSizes[encodedStream->size] = signalBits;
			encodedStream->elementSizes[encodedStream->size + 1] = signalBits;
			encodedStream->size += 2;
#endif
			for(int i = 0; i < subSegmentCount; i ++){
				if(i != oldMinPosition){
					transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 1));
					transmittedData[dataCount] = subMin[i];
					dataCount ++;
					////
					*bandwidth += 1 * dataUnitBits; // Extra max, min data points.
					////
#ifdef COMPRESSION
					encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 1));
					encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 1));
					encodedStream->stream[encodedStream->size] = subMin[i];
					encodedStream->elementSizes[encodedStream->size] = dataUnitBits;
					encodedStream->size += 1;
#endif
				}
				if(i != oldMaxPosition){
					transmittedData = (double*)realloc(transmittedData, sizeof(double) * (dataCount + 1));
					transmittedData[dataCount] = subMax[i];
					dataCount ++;
					////
					*bandwidth += 1 * dataUnitBits; // Extra max, min data points.
					////
#ifdef COMPRESSION
					encodedStream->stream = (double*)realloc(encodedStream->stream, sizeof(double) * (encodedStream->size + 1));
					encodedStream->elementSizes = (int*)realloc(encodedStream->elementSizes, sizeof(int) * (encodedStream->size + 1));
					encodedStream->stream[encodedStream->size] = subMax[i];
					encodedStream->elementSizes[encodedStream->size] = dataUnitBits;
					encodedStream->size += 1;
#endif
				}
			}
		}
	}
	return transmittedData;
}

// Decodes a data stream.
// sensor: A course time series at certain resolution level.
// transmittedData: The data stream.
// subSegmentCount: The parameter R. How many segments at the next level a segment at this level is cut into.
Sensor* sensorDecodeSplitInfo(Sensor *sensor, double *transmittedData, int subSegmentCount){
	int totalSegmentCount = 0;
	for(int t = 0; t < sensor->segmentCount; t ++){
		totalSegmentCount += (sensor->segments[t].size / subSegmentCount > 0) ? subSegmentCount : sensor->segments[t].size;
	}
	Sensor *newSensor = sensorConstruct(sensor->index, totalSegmentCount);
	
	for(int t = 0, v = 0, d = 0; t < sensor->segmentCount; t ++){
		// General rule 1
		if(sensor->segments[t].size == 1){
			newSensor->segments[v] = sensor->segments[t];
			v ++;
		}
		// General rule 2
		else if(sensor->segments[t].size == 2){
			double signal = transmittedData[d];
			d ++;
			if(signal == 0){
				newSensor->segments[v].min = sensor->segments[t].min;
				newSensor->segments[v].max = sensor->segments[t].min;
				newSensor->segments[v].size = 1;
				v ++;
				newSensor->segments[v].min = sensor->segments[t].max;
				newSensor->segments[v].max = sensor->segments[t].max;
				newSensor->segments[v].size = 1;
				v ++;
			}
			else{
				newSensor->segments[v].min = sensor->segments[t].max;
				newSensor->segments[v].max = sensor->segments[t].max;
				newSensor->segments[v].size = 1;
				v ++;	
				newSensor->segments[v].min = sensor->segments[t].min;
				newSensor->segments[v].max = sensor->segments[t].min;
				newSensor->segments[v].size = 1;
				v ++;
			}
		}
		// General rule 3
		else if(sensor->segments[t].size > 2 && sensor->segments[t].size <= subSegmentCount){
			int oldMinPosition = (int)transmittedData[d];
			int oldMaxPosition = (int)transmittedData[d + 1];
			d += 2;
			for(int i = 0; i < sensor->segments[t].size; i ++){
				if(i == oldMinPosition){
					newSensor->segments[v].min = sensor->segments[t].min;
					newSensor->segments[v].max = sensor->segments[t].min;
				}
				else if(i == oldMaxPosition){	
					newSensor->segments[v].min = sensor->segments[t].max;
					newSensor->segments[v].max = sensor->segments[t].max;
				}
				else{
					newSensor->segments[v].min = transmittedData[d];
					newSensor->segments[v].max = transmittedData[d];
					d ++;
				}
				newSensor->segments[v].size = 1;
				v ++;
			}
		}
		// General rule 4
		else if(sensor->segments[t].size > subSegmentCount && sensor->segments[t].size < subSegmentCount * 2){
			int oldMinPosition = (int)transmittedData[d];
			int oldMaxPosition = (int)transmittedData[d + 1];
			d += 2;
			for(int i = 0; i < subSegmentCount - 1; i ++){
				if(i == oldMinPosition){
					newSensor->segments[v].min = sensor->segments[t].min;
					newSensor->segments[v].max = sensor->segments[t].min;
				}
				else if(i == oldMaxPosition){
					newSensor->segments[v].min = sensor->segments[t].max;
					newSensor->segments[v].max = sensor->segments[t].max;
				}
				else{
					newSensor->segments[v].min = transmittedData[d];
					newSensor->segments[v].max = transmittedData[d];
					d ++;
				}
				newSensor->segments[v].size = 1;
				v ++;
			}
			if(subSegmentCount - 1 == oldMinPosition){
				newSensor->segments[v].min = sensor->segments[t].min;
			}
			else{
				newSensor->segments[v].min = transmittedData[d];
				d ++;
			}
			if(subSegmentCount - 1 == oldMaxPosition){
				newSensor->segments[v].max = sensor->segments[t].max;
			}
			else{
				newSensor->segments[v].max = transmittedData[d];
				d ++;
			}
			newSensor->segments[v].size = sensor->segments[t].size - (subSegmentCount - 1);
			v ++;
		}
		// General rule 5
		else if(sensor->segments[t].size >= subSegmentCount * 2){
			int oldMinPosition = (int)transmittedData[d];
			int oldMaxPosition = (int)transmittedData[d + 1];
			d += 2;
			int subSize = sensor->segments[t].size / subSegmentCount;
			int remainingSize = sensor->segments[t].size - subSize * (subSegmentCount - 1);
			for(int i = 0; i < subSegmentCount; i ++){
				if(i == oldMinPosition){
					newSensor->segments[v].min = sensor->segments[t].min;
				}
				else{
					newSensor->segments[v].min = transmittedData[d];
					d ++;
				}
				if(i == oldMaxPosition){
					newSensor->segments[v].max = sensor->segments[t].max;
				}
				else{
					newSensor->segments[v].max = transmittedData[d];
					d ++;
				}
				if(i < subSegmentCount - 1){
					newSensor->segments[v].size = subSize;
				}
				else{
					newSensor->segments[v].size = remainingSize;
				}
				v ++;
			}
		}
	}
	
	return newSensor;
}

// Detects whether a coarse query becomes fine. That is, asks whether Q^{l} = Q.
bool sensorIsAccurate(Sensor *sensor){
	for(int t = 0; t < sensor->segmentCount; t ++){
		if(sensor->segments[t].size > 1){
			return false;
		}
	}
	return true;
}

// Repeats data points to construct a coarse query Q^{l}.
Sensor* sensorRepeatSegments(Sensor *sensor){
	int totalSegmentCount = 0;
	for(int t = 0; t < sensor->segmentCount; t ++){
		totalSegmentCount += sensor->segments[t].size;
	}
	Sensor *newSensor = sensorConstruct(sensor->index, totalSegmentCount);

	for(int t = 0, startIndex = 0; t < sensor->segmentCount; t ++){
		for(int i = 0; i < sensor->segments[t].size; i ++){
			newSensor->segments[startIndex + i].min = sensor->segments[t].min;
			newSensor->segments[startIndex + i].max = sensor->segments[t].max;
			newSensor->segments[startIndex + i].size = 1;
		}
		startIndex += sensor->segments[t].size;
	}	
	return newSensor;
}

Sensor* sensorConstruct(int index, int realSegmentCount){
	Sensor* sensor = (Sensor*)malloc(sizeof(Sensor));
	sensor->index = index;
	sensor->realSegmentCount = realSegmentCount;
	sensor->segmentCount = realSegmentCount;
	sensor->realSegments = (Segment*)malloc(sizeof(Segment) * realSegmentCount);
	sensor->segments = (Segment*)malloc(sizeof(Segment) * realSegmentCount);
	sensor->coarseSensor = NULL;
	return sensor;
}

void sensorDestruct(Sensor* sensor){
	if(sensor->coarseSensor != NULL){
		if(sensor->coarseSensor->realSegments == sensor->coarseSensor->segments){
			sensor->coarseSensor->segments = NULL;
		}
		free(sensor->coarseSensor->realSegments);
		free(sensor->coarseSensor->segments);
		free(sensor->coarseSensor);
		if(sensor == sensor->coarseSensor){
			sensor = NULL;
		}
	}
	if(sensor != NULL){
		if(sensor->realSegments == sensor->segments){
			sensor->segments = NULL;
		}
		free(sensor->realSegments);
		free(sensor->segments);
		free(sensor);
	}
}

// Copies time series data.
Sensor* sensorCopy(Sensor *sensor){
	Sensor *newSensor = sensorConstruct(sensor->index, sensor->realSegmentCount);
	newSensor->segmentCount = sensor->segmentCount;
	for(int t = 0; t < sensor->segmentCount; t ++){
		newSensor->segments[t] = sensor->segments[t];
	}
	return newSensor;
}

//================
// Site-related functions.

Site* siteConstruct(int index){
	Site *site = (Site*)malloc(sizeof(Site));
	site->index = index;
	site->sensorCount = 0;
	site->sensors = NULL;
	return site;
}

void siteDestruct(Site *site){
	free(site->sensors);
	free(site);
	site = NULL;
}

//================
// Dataset-related functions.

// Read time series dataset. The file format follows UCR Time Series Datasets.
void datasetReadData(Dataset *dataset, char *dataFileName){
	int lineLength = 10000000;
	char *line = (char*)malloc(sizeof(char) * lineLength);
	dataset->sensors = (Sensor**)malloc(sizeof(Sensor*));
	FILE *dataFile = fopen(dataFileName, "r");

	for(dataset->realSensorCount = 0; fgets(line, lineLength, dataFile) != NULL; dataset->realSensorCount ++){
		int charCount;
		double value;
		bool isLabel = true;
		int realSegmentCount = 0;
		for(char *ptr = line; sscanf(ptr, "%lf%n", &value, &charCount) == 1; ptr += charCount){	
			if(isLabel){
				isLabel = false;
				
			}
			else{
				realSegmentCount ++;
			}
		}
		int s = dataset->realSensorCount;
		dataset->sensors = (Sensor**)realloc(dataset->sensors, sizeof(Sensor*) * (s + 1));
		dataset->sensors[s] = sensorConstruct(s, realSegmentCount);
		
		isLabel = true;
		int t = 0;
		for(char *ptr = line; sscanf(ptr, "%lf%n", &value, &charCount) == 1; ptr += charCount){	
			if(isLabel){
				isLabel = false;
				dataset->sensors[s]->label = value;
			}
			else{
				dataset->sensors[s]->realSegments[t].min = value;
				dataset->sensors[s]->realSegments[t].max = value;
				dataset->sensors[s]->realSegments[t].size = 1;
				t ++;
			}
		}	

	}
	fclose(dataFile);
	free(line);
}

// Constructs a dataset for k nearest neighbors.
// Equally distributes sensors to sites.
// sensorCount: The parameter S. The number of sensors participating the kNN.
// siteCount: The parameter M. The number of sites participating the kNN.
// queryLength: The length of the query and candidates.
// queryIndex: Assigns which sensor (time series) as the query.
void datasetSplitSites(Dataset *dataset, int sensorCount, int siteCount, int queryLength, int queryIndex){
	if(dataset->sites != NULL){
		for(int m = 0; m < dataset->siteCount; m ++){
			siteDestruct(dataset->sites[m]);
		}
		free(dataset->sites);
	}
	dataset->sites = (Site**)malloc(sizeof(Site*) * siteCount);
	for(int m = 0; m < siteCount; m ++){
		dataset->sites[m] = siteConstruct(m);
		dataset->sites[m]->sensors = (Sensor**)malloc(sizeof(Sensor*));
	}

	dataset->query = dataset->sensors[queryIndex];
	dataset->query->segmentCount = queryLength;
	dataset->query->segments = (Segment*)realloc(dataset->query->segments, sizeof(Segment) * queryLength);
	for(int t = 0; t < queryLength; t ++){
		dataset->query->segments[t] = dataset->query->realSegments[t];
	}
	sensorNormalizeSegments(dataset->query);

	int siteSensorCount = (sensorCount + siteCount - 1) / siteCount;
	for(int s = 0, m = 0; s < sensorCount + 1; s ++){
		if(s == queryIndex){
			continue;
		}
		int i = dataset->sites[m]->sensorCount;
		dataset->sites[m]->sensors = (Sensor**)realloc(dataset->sites[m]->sensors, sizeof(Sensor*) * (i + 1));
		dataset->sites[m]->sensors[i] = dataset->sensors[s];
		dataset->sites[m]->sensors[i]->segmentCount = queryLength;
		dataset->sites[m]->sensors[i]->segments = (Segment*)realloc(dataset->sites[m]->sensors[i]->segments, sizeof(Segment) * queryLength);
		for(int t = 0; t < queryLength; t ++){
			dataset->sites[m]->sensors[i]->segments[t] = dataset->sites[m]->sensors[i]->realSegments[t];
		}
		sensorNormalizeSegments(dataset->sites[m]->sensors[i]);
		dataset->sites[m]->sensorCount ++;
		m ++;
		if(m % siteCount == 0){
			m = 0;
		}
	}
	
	dataset->sensorCount = sensorCount;
	dataset->siteCount = siteCount;
	dataset->queryLength = queryLength;
}

// Returns the sensor instance given its index.
Sensor* datasetFindSensor(Dataset *dataset, int index){
	for(int m = 0; m < dataset->siteCount; m ++){
		Site *site = dataset->sites[m];
		for(int s = 0; s < site->sensorCount; s ++){
			Sensor *sensor = site->sensors[s];
			if(sensor->index == index){
				return sensor;
			}
		}
	}
}

// Returns the site instance given its index.
Site* datasetFindSite(Dataset *dataset, int index){
	for(int m = 0; m < dataset->siteCount; m ++){
		Site *site = dataset->sites[m];
		if(site->index == index){
			return site;
		}
	}
}

// Arranges the order of sensors at random.
void datasetRandomSensorOrder(Dataset *dataset){
	for(int i = 0; i < dataset->realSensorCount * 2; i ++){
		int randIndex1 = rand() % dataset->realSensorCount;
		int randIndex2 = rand() % dataset->realSensorCount;
		Sensor *temp = dataset->sensors[randIndex1];
		dataset->sensors[randIndex1] = dataset->sensors[randIndex2];
		dataset->sensors[randIndex2] = temp;
	}
}

Dataset* datasetConstruct(char *dataFileName){
	Dataset *dataset = (Dataset*)malloc(sizeof(Dataset));
	datasetReadData(dataset, dataFileName);
	dataset->siteCount = 0;
	dataset->sensorCount = dataset->realSensorCount;
	dataset->queryLength = 0;
	dataset->sites = NULL;
	return dataset;
}

void datasetDestruct(Dataset *dataset){
	for(int s = 0; s < dataset->realSensorCount; s ++){
		sensorDestruct(dataset->sensors[s]);
	}
	free(dataset->sensors);
	for(int m = 0; m < dataset->siteCount; m ++){
		siteDestruct(dataset->sites[m]);
	}
	free(dataset->sites);
	
	free(dataset);
}

//================
// Model (approach) -related functions.

int randomIntRange(int a, int b){
	if(a > b){
		int temp = a;
		a = b;
		b = temp;
	}
	return a + rand() % (b - a + 1);
}

// Increments the bandwidth.
// isData: true => data points having B bits.
void modelBandwidthIncrement(Model *model, int number, bool isData){
	if(model != NULL){
		model->bandwidth += number * ((isData)? model->dataUnitBits: 1);
	}
}

void modelResetBandwidth(Model *model){
	model->bandwidth = 0;
}

void modelPrint(Model *model, char *resultFileName){
	FILE *resultFile = fopen(resultFileName, "w");
	
	fprintf(resultFile, "Server: %d\t%.0lf\n", model->dataset->query->index, model->dataset->query->label);
	fprintf(resultFile, "K: %d\n", model->K);
	fprintf(resultFile, "Sensor Count: %d\n", model->dataset->sensorCount);
	fprintf(resultFile, "Site Count: %d\n", model->dataset->siteCount);
	fprintf(resultFile, "Query Length: %d\n", model->dataset->queryLength);
	fprintf(resultFile, "Bandwidth: %d\n", model->bandwidth);
	for(int i = 0; i < model->K; i ++){
		Sensor *sensor = datasetFindSensor(model->dataset, model->knnList[i].sensorIndex);
		fprintf(resultFile, "%d\t%.4lf\t%.0lf\n", model->knnList[i].sensorIndex, model->knnList[i].dtw, sensor->label);
	}

	fclose(resultFile);
}

int compareKNNDataDTW(const void *a, const void *b){
	KNNData *p = (KNNData*)a;
	KNNData *q = (KNNData*)b;

	if(p->dtw < q->dtw){
		return -1;
	}
	if(p->dtw > q->dtw){
		return 1;
	}
	if(p->sensorIndex < q->sensorIndex){
		return -1;
	}
	if(p->sensorIndex > q->sensorIndex){
		return 1;
	}
	return 0;
}

int compareKNNDataLowerBound(const void *a, const void *b){
	KNNData *p = (KNNData*)a;
	KNNData *q = (KNNData*)b;

	if(p->lowerBound < q->lowerBound){
		return -1;
	}
	if(p->lowerBound > q->lowerBound){
		return 1;
	}
	if(p->sensorIndex < q->sensorIndex){
		return -1;
	}
	if(p->sensorIndex > q->sensorIndex){
		return 1;
	}
	return 0;
}

void KNNDataRandomOrder(KNNData *knnData, int length){
	for(int r = 0; r < length * 2; r ++){
		int i = randomIntRange(0, length - 1);
		int j = randomIntRange(0, length - 1);
		KNNData temp = knnData[i];
		knnData[i] = knnData[j];
		knnData[j] = temp;
	}
}

// Calls one of lower bounds for DTW.
// type: Which lower bound called.
// threshold: For early abandoning. Negative threshold means no early abandoning.
double computeDTWBound(Sensor *query, Sensor *sensor, DTWLowerBound type, double threshold){
	if(threshold < 0){
		threshold = DBL_MAX;
	}
	double bound;
	switch(type){
		case LB_MS:
		{
			DTW *dtw;
			dtw = DTWConstruct(query, sensor);
			bound = DTWStart(dtw, threshold);
			DTWDestruct(dtw);
		}
			break;
		case LB_KEOGH:
		{
			LBKeogh *lbkeogh;
			lbkeogh = LBKeoghConstruct(query, sensor);
			double boundType1 = LBKeoghStart(lbkeogh, threshold);
			LBKeoghDestruct(lbkeogh);
			lbkeogh = LBKeoghConstruct(sensor, query);
			double boundType2 = LBKeoghStart(lbkeogh, threshold);
			LBKeoghDestruct(lbkeogh);
			
			bound = (boundType1 > boundType2) ? boundType1 : boundType2;
		}
			break;
		case LB_KIM:
		{
			LBKim *lbkim;
			lbkim = LBKimConstruct(query, sensor);
			bound = LBKimStart(lbkim);
			LBKimDestruct(lbkim);
		}
			break;
	}
	return bound;
}

bool exceedThreshold(double value, double threshold){
	return value - threshold > 0;
}

double BASELINE_COMPRESSION_RATE;

// Baseline approach.
void naiveKnnStart(Model *model, int K){
	if(model->knnList != NULL){
		free(model->knnList);
	}
	model->knnList = (KNNData*)malloc(sizeof(KNNData));
	modelResetBandwidth(model);
	
	Dataset *dataset = model->dataset;
	Sensor *query = dataset->query;
	int candidateCount = 0;
	
#ifdef COMPRESSION
	FileStream *baselineStream = fileStreamConstruct();
	baselineStream->size = dataset->queryLength;
	baselineStream->stream = (double*)realloc(baselineStream->stream, sizeof(double) * baselineStream->size);
	baselineStream->elementSizes = (int*)realloc(baselineStream->elementSizes, sizeof(int) * baselineStream->size);
	for(int t = 0; t < dataset->queryLength; t ++){
		baselineStream->stream[t] = query->segments[t].max;
		baselineStream->elementSizes[t] = model->dataUnitBits;
	}
	int queryOriginalSize, queryCompressedSize;
	BASELINE_COMPRESSION_RATE += fileStreamTestCompressionRate(baselineStream, &queryOriginalSize, &queryCompressedSize);
	fileStreamDestruct(baselineStream);
	printf("Baseline compression rate: %lf\n", BASELINE_COMPRESSION_RATE);
#endif

	for(int m = 0; m < dataset->siteCount; m ++){
		Site *site = dataset->sites[m];
		////
		modelBandwidthIncrement(model, dataset->queryLength, true); // Exact query (Server -> Site)
		////
		for(int s = 0; s < site->sensorCount; s ++){
			Sensor *sensor = site->sensors[s];
			double dtwValue = computeDTWBound(query, sensor, LB_MS, -1);

			model->knnList = (KNNData*)realloc(model->knnList, sizeof(KNNData) * (candidateCount + 1));
			model->knnList[candidateCount].sensorIndex = sensor->index;
			model->knnList[candidateCount].dtw = dtwValue;
			candidateCount ++;
			////
			modelBandwidthIncrement(model, 2, true); // Sensor index, DTW value (Site -> Server)
			////
		}
		////
		modelBandwidthIncrement(model, 1, true); // Site index (Site -> Server)
		////
	}
	qsort(model->knnList, candidateCount, sizeof(KNNData), compareKNNDataDTW);
	model->K = candidateCount;
}

int INITIAL_SENSOR_COUNT;
int FIRST_LEVEL_PRUNED_SENSOR_COUNT;
int LB_KIM_PRUNED_SENSOR_COUNT;
int LB_KEOGH_PRUNED_SENSOR_COUNT;
int LB_MS_PRUNED_SENSOR_COUNT;
int UNPRUNED_SENSOR_COUNT;
int SITE_LEVEL_COUNT[11];
double FRAMEWORK_COMPRESSION_RATE;

// Our framework.
void frameworkKnnStart(Model *model, int K){
	modelResetBandwidth(model);

	Dataset *dataset = model->dataset;
	Sensor *query = dataset->query;
	if(model->knnList != NULL){
		free(model->knnList);
	}
	model->knnList = NULL;
	
	// First level of the query (coarsest query).
	//int suggestedSegmentCount = (int)round(sqrt(query->segmentCount / 2.0)); // Rule of thumb: #segments = round(sqrt(length / 2))
	int suggestedSegmentCount = model->firstLevelSegmentCount;
#ifdef EQUAL_SIZE
	// Equal-sized segmentation
	int segmentSize = (query->segmentCount + suggestedSegmentCount - 1) / suggestedSegmentCount;
	query->coarseSensor = sensorSegment(query, &segmentSize, true);
#else
	// Unequal-sized segmentation
	int *segmentDistribution = sensorAnalyzeSegmentBottomUp(query, suggestedSegmentCount);
	query->coarseSensor = sensorSegment(query, segmentDistribution, false);
	free(segmentDistribution);
#endif
	////
	for(int m = 0; m < dataset->siteCount; m ++){
		for(int i = 0; i < suggestedSegmentCount; i ++){
			modelBandwidthIncrement(model, 3, true); // Initial min, max, size of segments (Server -> Site)
		}
		modelBandwidthIncrement(model, 1, true); // Number of sub-segments between two near levels (Server -> Site)
	}
	////

	// ================
	// Level-1 lower bounds (LB_MS).
	Sensor *firstLayerQuery = sensorCopy(query->coarseSensor);
	Sensor *reconstructedQuery = sensorRepeatSegments(query->coarseSensor);
	int knnListCount = 0;

#ifdef COMPRESSION
	int queryTotalOriginalSize, queryTotalCompressedSize;
	FileStream *level1Stream = fileStreamConstruct();
	level1Stream->size = reconstructedQuery->segmentCount * 3;
	level1Stream->stream = (double*)realloc(level1Stream->stream, sizeof(double) * level1Stream->size);
	level1Stream->elementSizes = (int*)realloc(level1Stream->elementSizes, sizeof(int) * level1Stream->size);
	for(int t = 0, i = 0; t < reconstructedQuery->segmentCount; t ++, i += 3){
		level1Stream->stream[i] = reconstructedQuery->segments[t].min;
		level1Stream->stream[i + 1] = reconstructedQuery->segments[t].max;
		level1Stream->stream[i + 2] = reconstructedQuery->segments[t].size;
		level1Stream->elementSizes[i] = model->dataUnitBits;
		level1Stream->elementSizes[i + 1] = model->dataUnitBits;
		level1Stream->elementSizes[i + 2] = model->dataUnitBits;
	}
	int queryOriginalSize, queryCompressedSize;
	fileStreamTestCompressionRate(level1Stream, &queryOriginalSize, &queryCompressedSize);
	fileStreamDestruct(level1Stream);
	queryTotalOriginalSize += queryOriginalSize * dataset->siteCount;
	queryTotalCompressedSize += queryCompressedSize * dataset->siteCount;
#endif

	for(int m = 0; m < dataset->siteCount; m ++){
		Site *site = dataset->sites[m];
		site->candidateCount = site->sensorCount;

		for(int s = 0; s < site->sensorCount; s ++){
			Sensor *sensor = site->sensors[s];
			sensor->isCandidate = true;

#ifdef ORACLE_INITIALIZATION	
			double lowerBound = computeDTWBound(query, sensor, LB_MS, -1);
#else
			double lowerBound = computeDTWBound(reconstructedQuery, sensor, LB_MS, -1);
#endif
			
			model->knnList = (KNNData*)realloc(model->knnList, sizeof(KNNData) * (knnListCount + 1));
			model->knnList[knnListCount].sensorIndex = sensor->index;
			model->knnList[knnListCount].siteIndex = site->index;
			model->knnList[knnListCount].lowerBound = lowerBound;
			knnListCount ++;
			////
			modelBandwidthIncrement(model, 2, true); // Sensor index, lower bound value (Site -> Server)
			////
		}
		////
		modelBandwidthIncrement(model, 1, true); // Site index (Site -> Server)
		////
	}
	sensorDestruct(reconstructedQuery);
	
	// ================
	// Initialization
	// Computing exact DTWs from sites that have the top-K lower bounds
#ifdef RANDOM_INITIALIZATION
	KNNData randomSites[dataset->siteCount];
	for(int m = 0; m < dataset->siteCount; m ++){
		randomSites[m].siteIndex = dataset->sites[m]->index;
	}
	KNNDataRandomOrder(randomSites, dataset->siteCount);
#else
	qsort(model->knnList, knnListCount, sizeof(KNNData), compareKNNDataLowerBound);
#endif
	
	int dtwListCount = 0;
	KNNData *dtwList = (KNNData*)malloc(sizeof(KNNData));

#ifdef COMPRESSION
	int initializedSiteCount = 0;
#endif

	for(int i = 0; i < K; i ++){
#ifdef RANDOM_INITIALIZATION
		if(i >= dataset->siteCount){
			break;
		}
		Site *initialSite = datasetFindSite(dataset, randomSites[i].siteIndex);
#else
		Site *initialSite = datasetFindSite(dataset, model->knnList[i].siteIndex);
#endif
		if(initialSite->candidateCount > 0){
#ifdef COMPRESSION
			initializedSiteCount ++;
#endif
			printf("Initial site: %d\n", initialSite->index);
			////
			modelBandwidthIncrement(model, query->segmentCount, true); // Exact query (Server -> Site)
			////
			for(int s = 0; s < initialSite->sensorCount; s ++){
				Sensor *sensor = initialSite->sensors[s];
				double dtwValue = computeDTWBound(query, sensor, LB_MS, -1);
				
				dtwList = (KNNData*)realloc(dtwList, sizeof(KNNData) * (dtwListCount + 1));
				dtwList[dtwListCount].sensorIndex = sensor->index;
				dtwList[dtwListCount].siteIndex = initialSite->index;
				dtwList[dtwListCount].dtw = dtwValue;
				dtwList[dtwListCount].lowerBound = dtwValue;
				dtwListCount ++;

				////
				modelBandwidthIncrement(model, 2, true); // Senser index, DTW value (Site -> Server)
				////
			}
			////
			modelBandwidthIncrement(model, 1, true); // Site index (Site -> Server)
			////
			initialSite->candidateCount = 0;
			
			INITIAL_SENSOR_COUNT += initialSite->sensorCount;
			SITE_LEVEL_COUNT[0] ++; 
		}
	}	
	qsort(dtwList, dtwListCount, sizeof(KNNData), compareKNNDataDTW);

#ifdef COMPRESSION
	printf("Qeury length: %d\n", dataset->queryLength);
	FileStream *initializationStream = fileStreamConstruct();
	initializationStream->size = dataset->queryLength;
	initializationStream->stream = (double*)realloc(initializationStream->stream, sizeof(double) * initializationStream->size);
	initializationStream->elementSizes = (int*)realloc(initializationStream->elementSizes, sizeof(int) * initializationStream->size);
	for(int t = 0; t < dataset->queryLength; t ++){
		initializationStream->stream[t] = query->segments[t].max;
		initializationStream->elementSizes[t] = model->dataUnitBits;
	}
	fileStreamTestCompressionRate(initializationStream, &queryOriginalSize, &queryCompressedSize);
	fileStreamDestruct(initializationStream);
	queryTotalOriginalSize += queryOriginalSize * initializedSiteCount;
	queryTotalCompressedSize += queryCompressedSize * initializedSiteCount;
#endif

	double threshold = dtwList[K - 1].dtw;
	dtwListCount = K;
	printf("Intiial threshold: %lf\n", threshold);

	// ================
	// Prunes sites using level-1 lower bounds
	KNNData minLowerBound[dataset->siteCount];
	for(int m = 0; m < dataset->siteCount; m ++){
		minLowerBound[m].siteIndex = dataset->sites[m]->index;
		minLowerBound[m].lowerBound = -1;
	}
	for(int i = 0; i < knnListCount; i ++){
		int index;
		for(int m = 0; m < dataset->siteCount; m ++){
			if(minLowerBound[m].siteIndex == model->knnList[i].siteIndex){
				index = m;
				break;
			}
		}

		if(minLowerBound[index].lowerBound == -1 || minLowerBound[index].lowerBound > model->knnList[i].lowerBound){
			minLowerBound[index].lowerBound = model->knnList[i].lowerBound;
		}
		
	}
	KNNData iterationOrder[dataset->siteCount];
	for(int m = 0; m < dataset->siteCount; m ++){
		iterationOrder[m] = minLowerBound[m];
		if(exceedThreshold(minLowerBound[m].lowerBound, threshold)){
			Site *site = datasetFindSite(dataset, minLowerBound[m].siteIndex);
			site->candidateCount = 0;  // Site is pruned.

			FIRST_LEVEL_PRUNED_SENSOR_COUNT += site->sensorCount;
			SITE_LEVEL_COUNT[1] ++;
		}
	}

	qsort(iterationOrder, dataset->siteCount, sizeof(KNNData), compareKNNDataLowerBound);

	// ================
	// Searches kNN on other sites
	for(int m = 0; m < dataset->siteCount; m ++){
		Site *site = datasetFindSite(dataset, iterationOrder[m].siteIndex);
		if(site->candidateCount > 0){
			int knnListCount = 0;
			bool atFinalLevel = false;
			int resolutionLevel = 1;
			
			free(query->coarseSensor);
			query->coarseSensor = sensorCopy(firstLayerQuery);

			////
			modelBandwidthIncrement(model, 1, true); // Threshold for Pruning (Server -> Site)
			////
			
			do{
				if(resolutionLevel < 10){
					resolutionLevel ++;
				}
				int sampleBandwidth;
				if(!atFinalLevel){
#ifdef COMPRESSION
					FileStream *sampleStream = fileStreamConstruct();
					double *transmittedData = sensorEncodeSplitInfo(query->coarseSensor, query, model->subSegmentCount, model->dataUnitBits, &sampleBandwidth, sampleStream);
					fileStreamTestCompressionRate(sampleStream, &queryOriginalSize, &queryCompressedSize);
					fileStreamDestruct(sampleStream);
					queryTotalOriginalSize += queryOriginalSize;
					queryTotalCompressedSize += queryCompressedSize;
#else
					double *transmittedData = sensorEncodeSplitInfo(query->coarseSensor, query, model->subSegmentCount, model->dataUnitBits, &sampleBandwidth, NULL);
#endif
					Sensor *newCoarseQuery = sensorDecodeSplitInfo(query->coarseSensor, transmittedData, model->subSegmentCount);
					free(transmittedData);
					sensorDestruct(query->coarseSensor);
					query->coarseSensor = newCoarseQuery;
					////
					modelBandwidthIncrement(model, sampleBandwidth, false); // Sampled data points (Server -> Site)
					////
					
					if(sensorIsAccurate(query->coarseSensor)){
						atFinalLevel = true;
					}
				}

				Sensor *reconstructedQuery = sensorRepeatSegments(query->coarseSensor);
				knnListCount = 0;
				bool isNull = true;
				for(int s = 0; s < site->sensorCount; s ++){
					Sensor *sensor = site->sensors[s];
					if(sensor->isCandidate){
						// Cascade pruning structure
						double lbKim = computeDTWBound(reconstructedQuery, sensor, LB_KIM, -1);
						if(exceedThreshold(lbKim, threshold)){
							sensor->isCandidate = false;
							site->candidateCount --;
							LB_KIM_PRUNED_SENSOR_COUNT ++;
							continue;
						}
						
						double lbKeogh = computeDTWBound(reconstructedQuery, sensor, LB_KEOGH, threshold);
						if(exceedThreshold(lbKeogh, threshold)){
							sensor->isCandidate = false;
							site->candidateCount --;
							LB_KEOGH_PRUNED_SENSOR_COUNT ++;
							continue;
						}

						double lbMS = computeDTWBound(reconstructedQuery, sensor, LB_MS, threshold);
						if(exceedThreshold(lbMS, threshold)){
							sensor->isCandidate = false;
							site->candidateCount --;
							LB_MS_PRUNED_SENSOR_COUNT ++;
							continue;
						}

						if(atFinalLevel){
							model->knnList = (KNNData*)realloc(model->knnList, sizeof(KNNData) * (knnListCount + 1));
							model->knnList[knnListCount].sensorIndex = sensor->index;
							model->knnList[knnListCount].siteIndex = site->index;
							model->knnList[knnListCount].dtw = lbMS; // LB_MS == DTW for the query at final resolution level.
							knnListCount ++;

							UNPRUNED_SENSOR_COUNT ++;
							////
							modelBandwidthIncrement(model, 2, true); // Sensor index, exact DTW to server (Site -> Server)
							////
						}
					}
				}
				sensorDestruct(reconstructedQuery);

				////
				modelBandwidthIncrement(model, 1, true); // Site index (Site -> Server)
				modelBandwidthIncrement(model, 1, false); // Signal: Whether continue to send samples or not (Site -> Server)
				////
				if(site->candidateCount == 0){
					atFinalLevel = false;
					
					SITE_LEVEL_COUNT[resolutionLevel] ++;
					break;
				}

			}while(!atFinalLevel);

			if(atFinalLevel){
				for(int i = 0; i < knnListCount; i ++){
					dtwList = (KNNData*)realloc(dtwList, sizeof(KNNData) * (dtwListCount + 1));
					dtwList[dtwListCount] = model->knnList[i];
					dtwListCount ++;	
				}
				qsort(dtwList, dtwListCount, sizeof(KNNData), compareKNNDataDTW);
				dtwList = (KNNData*)realloc(dtwList, sizeof(KNNData) * K);
				threshold = dtwList[K - 1].dtw;
				dtwListCount = K;
				printf("Site %d updates threshold to %.6lf\n", site->index, threshold);
				
				SITE_LEVEL_COUNT[resolutionLevel] ++;
			}
		}
	}
	free(model->knnList);
	model->knnList = dtwList;
	model->K = K;
#ifdef COMPRESSION
	FRAMEWORK_COMPRESSION_RATE += (double)queryTotalCompressedSize / queryTotalOriginalSize;
	//printf("Query total compressed size: %d\n", queryTotalCompressedSize);
	//printf("Query total original size: %d\n", queryTotalOriginalSize);
	printf("Framework compression rate: %lf\n", FRAMEWORK_COMPRESSION_RATE);
#endif
}

Model* modelConstruct(Dataset *dataset, void (*knnStart)(Model*, int), int firstLevelSegmentCount, int subSegmentCount){
	Model *model = (Model*)malloc(sizeof(Model));
	model->K = 0;
	model->bandwidth = 0;
	model->dataUnitBits = 32;
	model->firstLevelSegmentCount = firstLevelSegmentCount;
	model->subSegmentCount = subSegmentCount;
	model->knnList = NULL;
	model->dataset = dataset;
	model->knnStart = knnStart;
	return model;
}

void modelDestruct(Model *model){
	free(model->knnList);
	free(model);
	model = NULL;
}

//================

// Runs naive approach and framework.
bool testCode(Dataset *dataset, int sensorCount, int siteCount, int queryLength, int queryIndex, int K, int firstLevelSegmentCount, int subSegmentCount, int *bandwidthNaive, int *bandwidthFramework){	
	printf("Query index: %d\n", queryIndex);
	datasetSplitSites(dataset, sensorCount, siteCount, queryLength, queryIndex);
	printf("Dataset has been set\n");
	
	Model *naiveModel = modelConstruct(dataset, naiveKnnStart, firstLevelSegmentCount, subSegmentCount);
	naiveModel->knnStart(naiveModel, K);
	printf("Baseline Model has finished\n");
	
	Model *framework = modelConstruct(dataset, frameworkKnnStart, firstLevelSegmentCount, subSegmentCount);
	framework->knnStart(framework, K);
	printf("Framework Model has finished\n");
	
	bool modelCorrect = true;
	int correctCaseCount = 0;
	*bandwidthNaive = naiveModel->bandwidth;
	printf("\n");
	
	printf("Baseline:");
	for(int i = 0; i < K; i ++){
		printf("\t%d", naiveModel->knnList[i].sensorIndex);
	}
	printf("\n");
	printf("Framework:");
	for(int i = 0; i < K; i ++){
		printf("\t%d", framework->knnList[i].sensorIndex);
	}
	printf("\n");

	// Comparison of kNN between Naive and Framework
	for(int i = 0; i < K; i ++){
		bool found = false;
		for(int j = 0; j < K; j ++){
			if(naiveModel->knnList[j].sensorIndex == framework->knnList[i].sensorIndex){
				found = true;
				break;
			}
		}
		if(!found){
			modelCorrect = false;
		}
		else{
			correctCaseCount ++;
		}
	}
	printf("Ratio of correct cases: %d / %d\n", correctCaseCount, K);
	if(modelCorrect){
		printf("Framework correct!\n");
		*bandwidthFramework = framework->bandwidth;
		printf("\tBandwidth ratio: %d / %d = %lf\n", framework->bandwidth, naiveModel->bandwidth, (double)framework->bandwidth / naiveModel->bandwidth);
	}
	
	modelDestruct(naiveModel);
	modelDestruct(framework);

	return modelCorrect;
}

void showArgumentInstruction(){
	printf("-f: Dataset file name\n");
	printf("\n");

	printf("-k: The number of candidates to be retrieved (Default 10)\n");
	printf("\t1 <= k <= s\n");
	printf("\n");

	printf("-s: The number of sensors (candidates) (Default 500)\n");
	printf("\t1 <= s <= Total number of time series in a dataset - 1\n");
	printf("\tAble to test multiple values (e.g. 500,1000,1500) split by comma without space\n");
	printf("\n");

	printf("-m: The number of sites (Default 50)\n");
	printf("\t1 <= m <= s\n");
	printf("\tAble to test multiple values (e.g. 50,100,150) split by comma without space\n");
	printf("\n");

	printf("-t: Time series length (Default maximum length)\n");
	printf("\t1 <= t <= Maximum length provided by a dataset\n");
	printf("\n");

	printf("-q: Sensor (time series) index to be assigned as the query (Default -1)\n");
	printf("\t-2 <= q <= s\n");
	printf("\tq == -1: Let every time series be query in turn (for small data)\n");
	printf("\tq == -2: Randomly choose time series as query 100 times (for big data)\n");
	printf("\n");

	printf("-r: The parameter R. The number of segments cut from a segment at the previous level (Default 2)\n");
	printf("\t2 <= r\n");
	printf("\n");

	printf("-n: The parameter N. The number of segments at the first level (Default round(sqrt(t / 2)))\n");
	printf("\t1 <= n <= t\n");
}

int main(int argc, char *argv[]){
	char *datasetFileName = NULL;
	int K = 10;
	int sensorCountTestTimes = 1;
	int *sensorCounts = (int*)malloc(sizeof(int));
	sensorCounts[0] = 500;
	int siteCountTestTimes = 1;
	int *siteCounts = (int*)malloc(sizeof(int));
	siteCounts[0] = 50;
	int queryLength = -1;
	int queryIndex = -1;
	int subSegmentCount = 2;
	int firstLevelSegmentCount = -1;
	char *postfix = "";

	if(argc <= 1){
		showArgumentInstruction();
		return 0;
	}

	for(int a = 1; a < argc; a ++){
		if(argv[a][0] == '-'){
			switch(argv[a][1]){
				case 'f':
					datasetFileName = argv[a + 1];
					a ++;
					break;
				case 'k':
					sscanf(argv[a + 1], "%d", &K);
					a ++;
					break;
				case 's':
				{
					sensorCountTestTimes = 0;
					char *valueStr = strtok(argv[a + 1], ",");
					while(valueStr != NULL){
						sensorCounts = (int*)realloc(sensorCounts, sizeof(int) * (sensorCountTestTimes + 1));
						sscanf(valueStr, "%d", &sensorCounts[sensorCountTestTimes]);
						sensorCountTestTimes ++;
						valueStr = strtok(NULL, ",");
					}
					a ++;
				}
					break;
				case 'm':
				{
					siteCountTestTimes = 0;
					char *valueStr = strtok(argv[a + 1], ",");
					while(valueStr != NULL){
						siteCounts = (int*)realloc(siteCounts, sizeof(int) * (siteCountTestTimes + 1));
						sscanf(valueStr, "%d", &siteCounts[siteCountTestTimes]);
						siteCountTestTimes ++;
						valueStr = strtok(NULL, ",");
					}
					a ++;
				}
					break;
				case 't':
					sscanf(argv[a + 1], "%d", &queryLength);
					a ++;
					break;
				case 'q':
					sscanf(argv[a + 1], "%d", &queryIndex);
					a ++;
					break;
				case 'r':
					sscanf(argv[a + 1], "%d", &subSegmentCount);
					a ++;
					break;
				case 'n':
					sscanf(argv[a + 1], "%d", &firstLevelSegmentCount);
					a ++;
					break;
				case 'p':
					postfix = argv[a + 1];
					a ++;
					break;
			}
		}
		else{
			printf("Incorrect format of input arguments\n");
			return 0;
		}
	}

	if(datasetFileName == NULL){
		printf("No dataset file\n");
		return 0;
	}
	
	Dataset *dataset = datasetConstruct(datasetFileName);

	for(int s = 0; s < sensorCountTestTimes; s ++){
		for(int m = 0; m < siteCountTestTimes; m ++){
			if(siteCounts[m] <= 0 || siteCounts[m] > sensorCounts[s]){
				printf("Violate 1 <= Number of sites[%d] <= Number of sensors[%d]\n", m + 1, s + 1);
				return 0;
			}
		}
		if(K <= 0 || K > sensorCounts[s]){
			printf("Violate 1 <= k <= Number of sensors[%d]\n", s + 1);
			return 0;
		}
		if(queryIndex < -2 || queryIndex > sensorCounts[s]){
			printf("Violate -2 <= Sensor index <= Number of sensors[%d]\n", s + 1);
			return 0;
		}
		if(sensorCounts[s] <= 0){
			printf("Violate 1 <= Number of sensors[%d]\n", s + 1);
			return 0;
		}
	}

	int minStreamLength = -1;
	for(int s = 0; s < dataset->realSensorCount; s ++){
		if(minStreamLength == -1 || minStreamLength > dataset->sensors[s]->realSegmentCount){
			minStreamLength = dataset->sensors[s]->realSegmentCount;
		}
	}
	if(queryLength < 0){
		queryLength = minStreamLength;
	}
	if(queryLength <= 0 || queryLength > minStreamLength){
		printf("1 <= t <= Maximum length provided by a dataset\n");
		return 0;
	}

	if(firstLevelSegmentCount < 0){
		firstLevelSegmentCount = (int)round(sqrt(queryLength / 2.0));
	}
	if(firstLevelSegmentCount < 1 || firstLevelSegmentCount > queryLength){
		printf("1 <= Number of segments at the first level <= Time series length\n");
		return 0;
	}

	if(subSegmentCount < 2){
		printf("2 <= Number of segments cut from a segments at the previous level\n");
		return 0;
	}
	
	if(queryIndex < 0){
		int inFileNameLength = strlen(datasetFileName);

		char originalFileName[1000];
		strcpy(originalFileName, datasetFileName);
		char *originalDotPointer = strrchr(originalFileName, '.');
		if(originalDotPointer != NULL){
			originalDotPointer[0] = 0;
		}
		sprintf(originalFileName, "%s_original_file.out", originalFileName);
		ORIGINAL_FILE_NAME = originalFileName;
		
		char compressedFileName[1000];
		strcpy(compressedFileName, datasetFileName);
		char *compressedDotPointer = strrchr(compressedFileName, '.');
		if(compressedDotPointer != NULL){
			compressedDotPointer[0] = 0;
		}
		sprintf(compressedFileName, "%s_compressed_file.zip", compressedFileName);
		COMPRESSED_FILE_NAME = compressedFileName;

		double lastBandwidthRatioMean = 1;
		int bandwidth[2];
		int iterationTimes = 100;

		double bandwidthRatioMean[sensorCountTestTimes][siteCountTestTimes];
		double bandwidthRatioSTD[sensorCountTestTimes][siteCountTestTimes];
		double initialSensorCountMean[sensorCountTestTimes][siteCountTestTimes];
		double firstLevelPrunedSensorCountMean[sensorCountTestTimes][siteCountTestTimes];
		double lbKimPrunedSensorCountMean[sensorCountTestTimes][siteCountTestTimes];
		double lbKeoghPrunedSensorCountMean[sensorCountTestTimes][siteCountTestTimes];
		double lbMSPrunedSensorCountMean[sensorCountTestTimes][siteCountTestTimes];
		double unprunedSensorCountMean[sensorCountTestTimes][siteCountTestTimes];
		int resolutionLevelCount = 11;
		double siteLevelCountMean[sensorCountTestTimes][siteCountTestTimes][resolutionLevelCount];
		double baselineCompressionRate[sensorCountTestTimes][siteCountTestTimes];
		double frameworkCompressionRate[sensorCountTestTimes][siteCountTestTimes];
		
		for(int s = 0; s < sensorCountTestTimes; s ++){
			for(int m = 0; m < siteCountTestTimes; m ++){
				int siteCount = siteCounts[m];
				int sensorCount = sensorCounts[s];
				int siteSensorCount = sensorCount / siteCount;

				if(sensorCount > dataset->realSensorCount - 1 || K > sensorCount || siteSensorCount == 0){
					bandwidthRatioMean[s][m] = -1;
					bandwidthRatioSTD[s][m] = -1;
					initialSensorCountMean[s][m] = -1;
					firstLevelPrunedSensorCountMean[s][m] = -1;
					lbKimPrunedSensorCountMean[s][m] = -1;
					lbKeoghPrunedSensorCountMean[s][m] = -1;
					lbMSPrunedSensorCountMean[s][m] = -1;
					unprunedSensorCountMean[s][m] = -1;
					for(int p = 0; p < resolutionLevelCount; p ++){
						siteLevelCountMean[s][m][p] = -1;
					}
					baselineCompressionRate[s][m] = -1;
					frameworkCompressionRate[s][m] = -1;
					printf("The parameter combination is invalid for the dataset\n");
					continue;
				}
				
				bandwidthRatioMean[s][m] = 0;
				bandwidthRatioSTD[s][m] = 0;
				initialSensorCountMean[s][m] = 0;
				firstLevelPrunedSensorCountMean[s][m] = 0;
				lbKimPrunedSensorCountMean[s][m] = 0;
				lbKeoghPrunedSensorCountMean[s][m] = 0;
				lbMSPrunedSensorCountMean[s][m] = 0;
				unprunedSensorCountMean[s][m] = 0;
				for(int p = 0; p < resolutionLevelCount; p ++){
					siteLevelCountMean[s][m][p] = 0;
				}
				baselineCompressionRate[s][m] = 0;
				frameworkCompressionRate[s][m] = 0;

				srand(time(NULL) * lastBandwidthRatioMean);
				datasetRandomSensorOrder(dataset);

				if(queryIndex == -1){
					iterationTimes = sensorCount + 1;
				}

				for(int i = 0; i < iterationTimes; i ++){
					int thisQueryIndex;
					if(queryIndex == -1){
						thisQueryIndex = i;
					}
					else{
						thisQueryIndex = randomIntRange(0, sensorCount);
					}

					bool caseCorrect;
					
					INITIAL_SENSOR_COUNT = 0;
					FIRST_LEVEL_PRUNED_SENSOR_COUNT = 0;
					LB_KIM_PRUNED_SENSOR_COUNT = 0;
					LB_KEOGH_PRUNED_SENSOR_COUNT = 0;
					LB_MS_PRUNED_SENSOR_COUNT = 0;
					UNPRUNED_SENSOR_COUNT = 0;
					for(int p = 0; p < resolutionLevelCount; p ++){
						SITE_LEVEL_COUNT[p] = 0;
					}
					BASELINE_COMPRESSION_RATE = 0;
					FRAMEWORK_COMPRESSION_RATE = 0;
					
					caseCorrect = testCode(dataset, sensorCount, siteCount, minStreamLength, thisQueryIndex, K, firstLevelSegmentCount, subSegmentCount, &bandwidth[0], &bandwidth[1]);
					
					double bandwidthRatio = (double)bandwidth[1] / bandwidth[0];
					bandwidthRatioMean[s][m] += bandwidthRatio;
					bandwidthRatioSTD[s][m] += bandwidthRatio * bandwidthRatio;
					
					initialSensorCountMean[s][m] += INITIAL_SENSOR_COUNT;
					firstLevelPrunedSensorCountMean[s][m] += FIRST_LEVEL_PRUNED_SENSOR_COUNT;
					lbKimPrunedSensorCountMean[s][m] += LB_KIM_PRUNED_SENSOR_COUNT;
					lbKeoghPrunedSensorCountMean[s][m] += LB_KEOGH_PRUNED_SENSOR_COUNT;
					lbMSPrunedSensorCountMean[s][m] += LB_MS_PRUNED_SENSOR_COUNT;
					unprunedSensorCountMean[s][m] += UNPRUNED_SENSOR_COUNT;
					for(int p = 0; p < resolutionLevelCount; p ++){
						siteLevelCountMean[s][m][p] += SITE_LEVEL_COUNT[p];
					}
					baselineCompressionRate[s][m] += BASELINE_COMPRESSION_RATE;
					frameworkCompressionRate[s][m] += FRAMEWORK_COMPRESSION_RATE;

					printf("========\n");
				}
				bandwidthRatioMean[s][m] /= iterationTimes;
				bandwidthRatioSTD[s][m] = sqrt(bandwidthRatioSTD[s][m] / iterationTimes - bandwidthRatioMean[s][m] * bandwidthRatioMean[s][m]);
				initialSensorCountMean[s][m] /= iterationTimes;
				firstLevelPrunedSensorCountMean[s][m] /= iterationTimes;
				lbKimPrunedSensorCountMean[s][m] /= iterationTimes;
				lbKeoghPrunedSensorCountMean[s][m] /= iterationTimes;
				lbMSPrunedSensorCountMean[s][m] /= iterationTimes;
				unprunedSensorCountMean[s][m] /= iterationTimes;
				for(int p = 0; p < resolutionLevelCount; p ++){
					siteLevelCountMean[s][m][p] /= iterationTimes;
				}
				baselineCompressionRate[s][m] /= iterationTimes;
				frameworkCompressionRate[s][m] /= iterationTimes;

				lastBandwidthRatioMean = bandwidthRatioMean[s][m];
			}
		}

		printf("======== Prepare to write file ========\n");
		char outFileName[1000];
		strcpy(outFileName, datasetFileName);
		char *outDotPointer = strrchr(outFileName, '.');
		if(outDotPointer != NULL){
			outDotPointer[0] = 0;
		}
		sprintf(outFileName, "%s.out", outFileName);

		FILE *outFile = fopen(outFileName, "w");
		printf("The file \"%s\" is written\n", outFileName);
		int withoutLevelCount = 10;
		int itemCount = withoutLevelCount + resolutionLevelCount;
		for(int b = 0; b < itemCount; b ++){
			for(int m = -1; m < siteCountTestTimes; m ++){
				if(m == -1){
					fprintf(outFile, "%d", b + 1);
				}
				else{
					fprintf(outFile, "\t%d", siteCounts[m]);
				}
			}
			fprintf(outFile, "\n");
			for(int s = 0; s < sensorCountTestTimes; s ++){
				for(int m = -1; m < siteCountTestTimes; m ++){
					if(m == -1){
						fprintf(outFile, "%d", sensorCounts[s]);
					}
					else{
						switch(b){
							case 0:
								if(bandwidthRatioMean[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", bandwidthRatioMean[s][m]);
								}
								break;
							case 1:
								if(bandwidthRatioSTD[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", bandwidthRatioSTD[s][m]);
								}
								break;
							case 2:
								if(initialSensorCountMean[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", initialSensorCountMean[s][m]);
								}
								break;
							case 3:
								if(firstLevelPrunedSensorCountMean[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", firstLevelPrunedSensorCountMean[s][m]);
								}
								break;
							case 4:
								if(lbKimPrunedSensorCountMean[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", lbKimPrunedSensorCountMean[s][m]);
								}
								break;
							case 5:
								if(lbKeoghPrunedSensorCountMean[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", lbKeoghPrunedSensorCountMean[s][m]);
								}
								break;
							case 6:
								if(lbMSPrunedSensorCountMean[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", lbMSPrunedSensorCountMean[s][m]);
								}
								break;
							case 7:
								if(unprunedSensorCountMean[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", unprunedSensorCountMean[s][m]);
								}
								break;
							case 8:
								if(baselineCompressionRate[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", baselineCompressionRate[s][m]);
								}	
								break;
							case 9:
								if(frameworkCompressionRate[s][m] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", frameworkCompressionRate[s][m]);
								}	
								break;
							default:
								if(siteLevelCountMean[s][m][b - withoutLevelCount] == -1){
									fprintf(outFile, "\tN/A");
								}
								else{
									fprintf(outFile, "\t%.4lf", siteLevelCountMean[s][m][b - withoutLevelCount]);
								}
								break;
						}
					}
				}
				fprintf(outFile, "\n");
			}
			fprintf(outFile, "\n");
		}
		for(int b = 0; b < itemCount; b ++){
			fprintf(outFile, "%d", b + 1);
			for(int s = 0; s < sensorCountTestTimes; s ++){
				for(int m = 0; m < siteCountTestTimes; m ++){
					switch(b){
						case 0:
							if(bandwidthRatioMean[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", bandwidthRatioMean[s][m]);
							}
							break;
						case 1:
							if(bandwidthRatioSTD[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", bandwidthRatioSTD[s][m]);
							}
							break;
						case 2:
							if(initialSensorCountMean[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", initialSensorCountMean[s][m]);
							}
							break;
						case 3:
							if(firstLevelPrunedSensorCountMean[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", firstLevelPrunedSensorCountMean[s][m]);
							}
							break;
						case 4:
							if(lbKimPrunedSensorCountMean[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", lbKimPrunedSensorCountMean[s][m]);
							}
							break;
						case 5:
							if(lbKeoghPrunedSensorCountMean[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", lbKeoghPrunedSensorCountMean[s][m]);
							}
							break;
						case 6:
							if(lbMSPrunedSensorCountMean[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", lbMSPrunedSensorCountMean[s][m]);
							}
							break;
						case 7:
							if(unprunedSensorCountMean[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", unprunedSensorCountMean[s][m]);
							}
							break;
						case 8:
							if(baselineCompressionRate[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", baselineCompressionRate[s][m]);
							}
							break;
						case 9:
							if(frameworkCompressionRate[s][m] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", frameworkCompressionRate[s][m]);
							}
							break;
						default:
							if(siteLevelCountMean[s][m][b - withoutLevelCount] == -1){
								fprintf(outFile, "\tN/A");
							}
							else{
								fprintf(outFile, "\t%.4lf", siteLevelCountMean[s][m][b - withoutLevelCount]);
							}
							break;
					}
				}
			}
			fprintf(outFile, "\n");
		}
		fclose(outFile);
	}
	else{
		if(sensorCounts[0] > dataset->realSensorCount - 1){
			printf("Violate number of sensors <= Total number of time series in a dataset - 1\n");
		}
		int bandwidth[2];
		testCode(dataset, sensorCounts[0], siteCounts[0], queryLength, queryIndex, K, firstLevelSegmentCount, subSegmentCount, &bandwidth[0], &bandwidth[1]);
	}

	datasetDestruct(dataset);
	free(sensorCounts);
	free(siteCounts);
	return 0;
}
