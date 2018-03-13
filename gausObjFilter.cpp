#include <iostream>
#include 'Vec.h'
#include <vector>

// Helper
vector<vector<int>> createGaussian(int size) {

	// Mid point (max at v[max][max])
    int max = size/2;
  
  	//.initialize kernel
    std::vector<std::vector<int>> v(size,std::vector<int>(size));
	for (int i = 0; i < size; ++i)
        for(int j = 0; j < size; ++j)
            v[i][j] = 0;
    
    // set first row (Pascal's Triangle)
    int x = 1;
    for (int i = 0; i <= max; ++i) {
        v[0][i] = x;
        x = x * ((size-1) - i) / (i + 1);
    }
    
    // fill gaussian  
    for (int j = 1; j <= max; ++j) 
        for (int i = j; i <= max; ++i)
            v[j][i] = v[0][i] * v[0][j];
            
    // flip diagonally in first quadrant
    for (int j = 0; j <= max; ++j) 
        for (int i = j; i <= max; ++i)             
            v[i][j] = v[j][i];
    
    // flip accross vertical line in upper half
    for (int j = 0; j <= max; ++j) 
        for (int i = 0; i < max; ++i) 
            v[j][size-i-1] = v[j][i];
    
    // flip accross horizontal line in lower half
    for (int i = 0; i < size; ++i) 
        for (int j = 0; j < max; ++j) 
            v[size-j-1][i] = v[j][i];
    
    // DEBUG
    /*
    for (int j = 0; j < size; ++j) {
        for (int i = 0; i < size; ++i) {
            std::cout << v[j][i] << ' ';       
        }        
        std::cout << std::endl;
    }
    */
    return v;
}

Vec gausObjFilter(std::vector<Vec> c, std::vector<int> o, int h, int w, int size) {

	// size of filter should be odd
	if (size % 2 != 0) {
		++size;
	}

	// create the gaussian filter
	std::vector<std::vector<int>> kernel = createGaussian(size);

	// total size of image
	int pix = h * w;
	
	// initialize return vector
	std::vector<Vec> g(pix);
	
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
						pixelSum = pixelSum + c[pi] * kernel[v][u];
						pixelNum += kernel[v][u];
					}
				}			
			}
		} // filterPix has all contributing pixel values
	
		// Compute the median value
		g[i] = pixelSum / (double)pixelNum;
	}
	
	// Return the vector of median values
	return g;		

}