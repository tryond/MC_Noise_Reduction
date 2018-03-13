#include <vector>
#include 'Vec.h'

std::vector<Vec> boxObjFilter(std::vector<Vec> c, std::vector<int> o, int h, int w, int size) {

	// size of filter should be odd
	if (size % 2 != 0) {
		++size;
	}

	// total size of image
	int pix = h * w;
	
	// initialize return vector
	std::vector<Vec> b(pix);
	
	// corner move
	int move = (size-1)/2;
	
	// corner index
	int ci, cx, cy;
	
	// point indices
	int x, y;
		
	// potential x, y, and i
	int px, py, pi;
		
	// filter pixel sum
	Vec pixelSum;
	int pixelNum;
		
	// For each pixel in the source
	for (int i = 0; i < pix; ++i) {
		
		// Set Vec to center pixel
		pixelSum = c[i];
		pixelNum = 1;
		
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
						pixelSum = pixelSum + c[pi];
						++pixelNum;
					}
				}			
			}
		} // filterPix has all contributing pixel values
	
		// Compute the median value
		b[i] = pixelSum / (double)pixelNum;
	}
	
	// Return the vector of median values
	return b;			
}