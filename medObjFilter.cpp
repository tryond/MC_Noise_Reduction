// Add magnitude to Vec
// Add comparator to Vec

#include <vector>
#include <cmath>
#include <algorithm>
#include 'Vec.h'

// ********************************
// Consider implementing median of medians here

// Helper function
Vec medianVec(std::vector<Vec> v) {

	sort(v.begin(), v.end());
	
	if (v.size() % 2 == 0)
		return v[v.size()/2].avg(v[v.size/2-1]);

	return v[v.size()/2];
}

std::vector<Vec> medObjFilter(std::vector<Vec> c, std::vector<int> o, int h, int w, int size) {

	// size of filter should be odd
	if (size % 2 != 0) {
		++size;
	}

	// total size of image
	int pix = h * w;
	
	// initialize return vector
	std::vector<Vec> m(pix);
	
	// corner move
	int move = (size-1)/2;
	
	// corner index
	int ci, cx, cy;
	
	// point indices
	int x, y;
		
	// potential x, y, and i
	int px, py, pi;
		
	// For each pixel in the source
	for (int i = 0; i < pix; ++i) {
		
		// Set Vec to center pixel
		Vec filterPix(c[i]);
		
		// Find x and y values
		x = i % w;
		y = i / w;
		
		// Find corner indices
		cx = x - move;
		cy = y - move;
	
		for (int v = 0; v < size; ++v) {
			for (int u = 0; u < size; ++u) {
			
				// Get potential pixel's x and y coords
				px = u + cx;
				py = v + cy;
				
				// Check that index is valid
				if (px >= 0 && px < w && py >= 0 && py < h) {
					
					// Compute potential pixel's index
					pi = (h - cy - 1) * w + cx;
					
					// Check that object id's are the same
					if (o[pi] == o[i]) {
						filterPix.push_back(c[pi]);
					}
				}			
			}
		} // filterPix has all contributing pixel values
	
		// Compute the median value
		m[i] = medianVec(filterPix);
	}
	
	// Return the vector of median values
	return m;			
}