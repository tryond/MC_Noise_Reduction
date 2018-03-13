#include "DofFilter.h"

std::vector<Vec> dofGray(std::vector<Vec> c, std::vector<int> o, std::vector<int> d, int h, int w, int levels, double focLevel) {

    if (levels <= 1)
        return c;

    // if focus too large or left undefined, set to middle
    if (focLevel < 0 || focLevel >= levels) {
        focLevel = floor((double)levels/2.0);
    }

    // DEBUG
    std::cout << "dofGray: " << levels << " / " << focLevel << std::endl;

    // total size of image
    int pix = h * w;
    
    // initialize return vector
    std::vector<Vec> b(pix);

    // Filter size for each pixel
    std::vector<int> s(pix);
      
    // corner move
    int move;
    
    // corner index
    int ci, cx, cy;
    
    // point indices
    int x, y;
        
    // potential x, y, and i
    int px, py, pi;
        
    // filter pixel sum
    Vec pixelSum;
    int pixelNum;
    
    // filter size (per pixel)
    int size;

    // Get filter size for each pixel
    double dmin = std::numeric_limits<double>::max(); 
    double dmax = std::numeric_limits<double>::min();
    
    // if focal distance not set, set to midpoint of scene
    double total = 0;
    for (int i = 0; i < pix; ++i) {
        if (d[i] < dmin) { dmin = d[i]; }
        if (d[i] > dmax) { dmax = d[i]; }
        total += d[i];
    }
    

    // Set boundaries according to levels
    double binSize = (dmax-dmin)/(double)levels;

    // How many bins on either side of focLevel
    // int focLevel = focLevel; // floor((focLevel-dmin)/binSize);

    std::vector<int> levelFilter(levels);

    /*
    std::cout << "dmin: " << dmin << std::endl;
    std::cout << "dmax: " << dmax << std::endl;
    std::cout << "focLevel: " << focLevel << std::endl;
    std::cout << "binSize: " << binSize << std::endl;
    std::cout << "focLevel: " << focLevel << std::endl;
    */

    // fill lower
    int filterSize = -1;
    for (int i = focLevel; i >= 0; --i) {
        filterSize += 2;
        levelFilter[i] = filterSize;
    }
    // fill upper
    filterSize = 1;
    for (int i = focLevel+1; i < levels; ++i) {
        filterSize += 2;
        levelFilter[i] = filterSize;
    }

    for (int i = 0; i < levels; ++i)
        std::cout << levelFilter[i] << ' ';
    std::cout << std::endl;

    // fill pixel filter sizes
    int binNo;
    for (int i = 0; i < pix; ++i) {
        binNo = floor((d[i]-dmin)/binSize);
        s[i] = levelFilter[binNo];
    }

    int unique_filters = levels - focLevel > focLevel ? levels - focLevel : focLevel;

    Vec green(0.0,1.0,0.0);
    Vec red(1.0,0.0,0.0);
    Vec blue(0.0,0.0,0.25);

    // For each pixel in the source
    for (int i = 0; i < pix; ++i) {
        double value = (double)(unique_filters - (((s[i] + 1.0) / 2.0) - 1.0)) * (1.0/(double)(unique_filters));
        Vec newVal = (green*value) + (red*(1.0-value)) + blue;
        b[i] = newVal;
    }
    
    return b;   // Return the vector of median values         
}

std::vector<Vec> dofMean(std::vector<Vec> c, std::vector<int> o, std::vector<int> d, int h, int w, int levels, double focLevel) {

    if (levels <= 1)
        return c;

    // if focus too large or left undefined, set to middle
    if (focLevel < 0 || focLevel >= levels) {
        focLevel = floor((double)levels/2.0);
    }

    // DEBUG
    std::cout << "dofMean: " << levels << " / " << focLevel << std::endl;

    // total size of image
    int pix = h * w;
    
    // initialize return vector
    std::vector<Vec> b(pix);

    // Filter size for each pixel
    std::vector<int> s(pix);
      
    // corner move
    int move;
    
    // corner index
    int ci, cx, cy;
    
    // point indices
    int x, y;
        
    // potential x, y, and i
    int px, py, pi;
        
    // filter pixel sum
    Vec pixelSum;
    int pixelNum;
    
    // filter size (per pixel)
    int size;

    // Get filter size for each pixel
    double dmin = std::numeric_limits<double>::max(); 
    double dmax = std::numeric_limits<double>::min();
    
    // if focal distance not set, set to midpoint of scene
    double total = 0;
    for (int i = 0; i < pix; ++i) {
        if (d[i] < dmin) { dmin = d[i]; }
        if (d[i] > dmax) { dmax = d[i]; }
        total += d[i];
    }

    // Set boundaries according to levels
    double binSize = (dmax-dmin)/(double)levels;

    // How many bins on either side of focLevel
    // int focLevel = focLevel; // floor((focLevel-dmin)/binSize);

    std::vector<int> levelFilter(levels);

    /*
    std::cout << "dmin: " << dmin << std::endl;
    std::cout << "dmax: " << dmax << std::endl;
    std::cout << "focLevel: " << focLevel << std::endl;
    std::cout << "binSize: " << binSize << std::endl;
    std::cout << "focLevel: " << focLevel << std::endl;
    */

    // fill lower
    int filterSize = -1;
    for (int i = focLevel; i >= 0; --i) {
        filterSize += 2;
        levelFilter[i] = filterSize;
    }
    // fill upper
    filterSize = 1;
    for (int i = focLevel+1; i < levels; ++i) {
        filterSize += 2;
        levelFilter[i] = filterSize;
    }

    for (int i = 0; i < levels; ++i)
        std::cout << levelFilter[i] << ' ';
    std::cout << std::endl;

    // fill pixel filter sizes
    int binNo;
    for (int i = 0; i < pix; ++i) {
        binNo = floor((d[i]-dmin)/binSize);
        s[i] = levelFilter[binNo];
    }

    // For each pixel in the source
    for (int i = 0; i < pix; ++i) {
        
        // Set Vec to center pixel
        pixelSum = c[i];
        pixelNum = 1;
        
        // Set size of kernel
        size = s[i];
        move = (size-1)/2; 

        // Find x and y values
        x = i % w;
        y = ( (h*w) - i - 1) / w;
        
        // Find corner indices
        cx = x - move;
        cy = y - move;
    
        if (size == 1) {
            b[i] = c[i];
            continue;
        }

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
                        pixelSum = pixelSum + c[pi];
                        ++pixelNum;
                    }

                }           
            }
        } // filterPix has all contributing pixel values
    
        // Compute the median value
        b[i] = pixelSum * (1.0/(double)pixelNum);
    }
    
    return b;   // Return the vector of median values         
}