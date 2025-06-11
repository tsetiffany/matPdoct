# PDOCT Post-Processing code

------ VERSION HISTORY ------
 
## VERSION 2.X ## 
v2.0 `main` branch push 05/27/2025
**MAJOR UPDATES**
-   2-frame averaging (previously 4-frame) to fix doubling/blurriness in DOPU composite images
-   `genDOPu_combinedBscans.m`: adjust uint8 range [0.3 1] (previously [0 1]) before converting to RGB for better DOPU colour range
-   LUT update `LUT_5050`

## VERSION 1.X ## 
v1.0.1 `OCTA` branch push
- **Initial push**
  - Handles 4CM OCTA, 500x500 with the correct reordering of B-scans (0101 0101 2323 2323...)
    
v1.0
* functional code, saving octv.mat and dopu.mat 
* Sockeye is running this version of the code
