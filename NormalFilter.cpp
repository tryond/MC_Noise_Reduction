#include "NormalFilter.h"

std::vector<Vec> normalMean(std::vector<Vec> c, std::vector<int> o, std::vector<Vec> n, int h, int w, int size) {

    if (size <= 1)
        return c;

    // size of filter should be odd
    if (size % 2 == 0) {
        ++size;
    }

    // DEBUG
    std::cout << "normalMean: " << size << std::endl;

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
    double pixelNum;
        
    // For each pixel in the source
    for (int i = 0; i < pix; ++i) {
        
        // Set Vec to center pixel
        pixelSum = c[i];
        pixelNum = 1;
        
        // Find x and y values
        x = i % w;
        y = ( (h*w) - i - 1) / w;
        
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
                    pi = (h - py - 1) * w + px;
                    
                    
                    // Check that objects match
                    if (o[pi] == o[i]) {
                        // How similar are the normals
                    	double scalar = n[i].dot(n[pi]);
                    	pixelSum = pixelSum + (c[pi] * scalar);
                    	pixelNum += scalar;
                    }
                    

                    // How similar are the normals
                    // double scalar = n[i].dot(n[pi]);
                    // pixelSum = pixelSum + (c[pi] * scalar);
                    // pixelNum += scalar;
                }           
            }
        } // filterPix has all contributing pixel values
    
        // Compute the median value
        b[i] = pixelSum * (1.0/(double)pixelNum);
    }
    
    return b;	// Return the vector of median values         
}


std::vector<std::vector<int>> normalGaussianFilter(int size) {

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
    
    return v;
}

std::vector<Vec> normalGaussian(std::vector<Vec> c, std::vector<int> o, std::vector<Vec> n, int h, int w, int size) {

    if (size <= 1)
        return c;

    // size of filter should be odd
    if (size % 2 == 0) {
        ++size;
    }

    // DEBUG
    std::cout << "normalGaussian: " << size << std::endl;

    // create the gaussian filter
    std::vector<std::vector<int>> kernel = normalGaussianFilter(size);

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
    double pixelNum;
    double scalar1;
    double scalar2;

    // For each pixel in the source
    for (int i = 0; i < pix; ++i) {
        
        // Set Vec to center pixel
        pixelSum = c[i];
        pixelNum = 1;
        
        // Find x and y values
        x = i % w;
        y = ( (h*w) - i - 1) / w;
        
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
                    pi = (h - py - 1) * w + px;

                    // Check that object id's are the same
                    if (o[pi] == o[i]) {
                        scalar1 = c[i].dot(c[pi]);
                        scalar2 = kernel[v][u] * scalar1;
                        pixelSum = pixelSum + c[pi] * scalar2;
                        pixelNum += scalar2;;
                    }
                }           
            }
        } // filterPix has all contributing pixel values
    
        // Compute the median value
        g[i] = pixelSum * (1.0/(double)pixelNum);
    }
    return g;       
}
