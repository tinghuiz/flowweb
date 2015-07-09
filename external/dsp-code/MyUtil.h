#pragma once

#include <string.h>
#include <limits>

using namespace std;

// function object to return indices of sorted array
template<class T> struct index_cmp_less {
	index_cmp_less(const T& arr):arr(arr) {}
	bool operator()(const size_t a, const size_t b) const
	{ return arr[a] < arr[b]; }
	const T& arr;
};

template<class T> struct index_cmp_greater {
	index_cmp_greater(const T& arr):arr(arr) {}
	bool operator()(const size_t a, const size_t b) const
	{ return arr[a] > arr[b]; }
	const T& arr;
};


//template <class T>
//void NMinS(T* src_im, T* dst_im, int sz, int width, int height)
//// Non-Minimum Suppression
//{
//	memcpy(dst_im, src_im, sizeof(T)*width*height);
//
//	// along the column
//	T* tmp1 = new T[width*height];
//	memcpy(tmp, src_im, sizeof(T)*width*height);
//	T limit_max = numeric_limits<int>::max();
//	for ( int x=0; x < width ; x++) {
//		for ( int y=0; y < height; y++) {
//			int idx1 = y +  x*height;
//			T val1 = src_im[idx1];
//			for ( int yy = y-sz ; yy <= y+sz ; yy++) {
//				if ( yy < 0 | yy >= height)
//					continue;
//
//				int idx2 = yy +  x*height;
//				T val2 = src_im[idx2];
//				if ( val1 > val2 ) { // if val1 is not minimum
//					dst_im[idx1] = limit_max;
//					tmp1[idx1] = val2;					
//				}
//			}
//		}
//	}
//
//	// along the row
//	for ( int y=0 ; y < width ; y++) {
//		for ( int x=0 ; x < height; x++) {
//			int idx1 = y + x*height;
//			T val1 = src_im[idx1];
//
//			for ( int xx = x-sz ; xx <= x+sz ; xx++) {
//				if ( xx < 0 | xx >= width ) 
//					continue;
//
//				int idx2 = y + xx*height;
//				T val2 = tmp1[idx2];
//				if ( val1 > val2 ) { 
//					dst_im[idx1] = limit_max;
//				}
//			}
//		}
//	}
//	delete [] tmp1;
//}
//
//template <class T>
//void NMaxS(T* src_im, T* dst_im, int sz, int width, int height)
//// Non-Maximum Suppression
//{
//	memcpy(dst_im, src_im, sizeof(T)*width*height);
//
//	// along the column
//	T* tmp1 = new T[width*height];
//	memcpy(tmp, src_im, sizeof(T)*width*height);
//	T limit_min = numeric_limits<int>::min();
//	for ( int x=0; x < width ; x++) {
//		for ( int y=0; y < height; y++) {
//			int idx1 = y + x*height;
//			T val1 = src_im[idx1];
//			for ( int yy = y-sz ; yy <= y+sz ; yy++) {
//
//				if ( yy < 0 | yy >= height)
//					continue;
//
//				int idx2 = yy + x*height;
//				T val2 = src_im[idx2];
//				if ( val1 < val2 ) { // if val1 is not maximum
//					dst_im[idx1] = limit_min;
//					tmp1[idx1] = val2;					
//				}
//			}
//		}
//	}
//
//	// along the row
//	for ( int y=0 ; y < width ; y++) {
//		for ( int x=0 ; x < height; x++) {
//			int idx1 = y + x*height;
//			T val1 = src_im[idx1];
//
//			for ( int xx = x-sz ; xx <= x+sz ; xx++) {
//				if ( xx < 0 | xx >= width ) 
//					continue;
//
//				int idx2 = y + xx*height;
//				T val2 = tmp1[idx2];
//				if ( val1 < val2 ) { 
//					dst_im[idx1] = limit_min;
//				}
//			}
//		}
//	}
//	delete [] tmp1;
//}
//				








