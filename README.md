<div align="center">
  <img src="pic.jpeg" alt="boaz"  width="1000" height="250"/>
  <h1>Ex0 - GNSS Raw Measurements 🌍📝🔍</h1>
  <p>Demo your location</p>
</div>

---

# 

## Overview

This assignment focuses on the basic principles of Global Navigation Satellite System (GNSS). WE implemented a naive positioning algorithm based on the Root Mean Square (RMS) of selected (i.e., weighted) pseudo-ranges. based on [John Mitchell's Website](https://www.johnsonmitchelld.com/)


## Objectives

1. Understand the basic concepts of GNSS, particularly pseudo-ranges.
2. Parse a dataset of GNSS raw measurements into a structured format.
3. Implement a positioning algorithm to compute positions based on GNSS data.
4. Convert ECEF coordinates to latitude, longitude, and altitude.
5. Integrate these tasks into a complete solution that processes raw GNSS measurements and outputs the computed path and positions in specified formats.
6. Perform testing on the provided dataset and additional data files.

## Tasks

1. **Familiarize with GNSS Concepts**:
   - Study the basic concepts of GNSS, particularly the notion of pseudo-ranges.
   - Review the provided presentation, "GNSS Raw Measurements", and the Android GnssLogger App.

2. **Data Parsing**:
   - Download the provided dataset.
   - Design a parsing tool to convert the log file into a CSV file with the following format:
     ```
     GPS time, SatPRN (ID), Sat.X, Sat.Y, Sat.Z, Pseudo-Range, CN0, Doppler (optional)
     ```
   - Note: The satellite positions should be in ECEF coordinates.
   - Refer to this [link](insert link) for a complete explanation and example Python code to aid in this task.

3. **Positioning Algorithm**:
   - Implement a positioning algorithm that, given a GPS time, computes the position in X, Y, Z coordinates.
   - Use an iterative numerical minimal RMS algorithm on a weighted set of SatPRNs.

4. **Coordinate Conversion**:
   - Implement a method to convert ECEF coordinates (X, Y, Z) to latitude, longitude, and altitude.
   - Refer to the [wiki page](insert link) for guidance.

5. **Integration and Output**:
   - Integrate the parsing, positioning algorithm, and coordinate conversion into a complete solution.
   - The solution:
     1. Receive a raw GNSS measurements log file.
     2. Output a KML file with the computed path (including time and animation).
     3. Output a CSV file with the following additional columns:
        ```
        Pos.X, Pos.Y, Pos.Z, Lat, Lon, Alt
        ```

6. **Testing and Documentation**:
   - Perform testing using the provided dataset and additional data files.

## Prerequisites

- Python 3.6+
- numpy
- pandas
- pyproj

Install the required packages using:
```bash
pip install numpy pandas pyproj
```
## How to Run

 **Setup**:
   - Clone the repository:
     ```bash
     git clone <repo_url>
     cd <repo_directory>
     ```
To change the data change data.txt
For solotion run:
 ```bash
     python3 solution.py
    python3 test.py
```

Feel free to reach out if you have any questions or need further assistance. Happy coding!





