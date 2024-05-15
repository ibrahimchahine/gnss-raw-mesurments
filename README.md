# Ex0 - GNSS Raw Measurements

## Overview

This assignment focuses on the basic principles of Global Navigation Satellite System (GNSS). WE implemented a naive positioning algorithm based on the Root Mean Square (RMS) of selected (i.e., weighted) pseudo-ranges.

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

## How to Run

1. **Setup**:
   - Clone the repository:
     ```bash
     git clone <repo_url>
     cd <repo_directory>
     ```
   - Install required dependencies:
     ```bash
     pip install -r requirements.txt
     ```

2. **Data Parsing**:
   - Place your GNSS raw measurements log file in the `data/` directory.
   - Run the parser script to generate the CSV file:
     ```bash
     python parse_log.py --input data/log_file.txt --output data/parsed_data.csv
     ```

3. **Positioning Algorithm**:
   - Run the positioning algorithm on the parsed CSV file:
     ```bash
     python compute_position.py --input data/parsed_data.csv --output data/positioned_data.csv
     ```

4. **Coordinate Conversion**:
   - Convert ECEF coordinates to latitude, longitude, and altitude:
     ```bash
     python convert_coordinates.py --input data/positioned_data.csv --output data/final_data.csv
     ```

5. **Generate Output Files**:
   - Generate the KML file and final CSV file:
     ```bash
     python generate_output.py --input data/final_data.csv --kml_output output/path.kml --csv_output output/final_positions.csv
     ```

6. **Testing**:
   - Run the test scripts to validate your implementation:
     ```bash
     python test_all.py
     ```
---

Feel free to reach out if you have any questions or need further assistance. Happy coding!





