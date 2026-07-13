// Fiji/ImageJ macro
// Batch rolling-ball background subtraction for files matching patternString
//
// Input:
// Any .tif / .tiff file whose filename contains patternString
//
// Output:
// <original_name>_bgsub.tif
// rolling_ball_log.csv

requires("1.53");

// ---------------- USER SETTINGS ----------------

// Only process files whose name contains this string.
// Example: "C561", "C488", "_Nuclei_1", "aligned.ome"
patternString = "C640";

rollingRadius = 100;        // pixels; set 0 to skip rolling-ball background subtraction
saveProcessedImages = true;

// Use Bio-Formats non-interactive importer.
// Recommended for .ome.tif / .ome.tiff files.
// Do NOT use for regular tif, otherwise Z channel will be T channel
useBioFormats = false;

// ------------------------------------------------

inputDir = getDirectory("Choose folder containing channel TIFFs");
outputDir = getDirectory("Choose output folder");

// ---------- Create processing log ----------
logPath = outputDir + "rolling_ball_log.csv";
File.saveString(
    "unique,inputFile,patternString,rollingRadius_px,backgroundSubtractionApplied,outputFile,status\n",
    logPath
);

files = getFileList(inputDir);

setBatchMode(true);

for (i = 0; i < files.length; i++) {

    fileName = files[i];

    // Skip folders.
    if (File.isDirectory(inputDir + fileName)) {
        continue;
    }

    // Only process tif/tiff files.
    if (!(endsWith(fileName, ".tif") || endsWith(fileName, ".tiff"))) {
        continue;
    }

    // Only process files containing patternString.
    if (indexOf(fileName, patternString) < 0) {
        continue;
    }

    pathA = inputDir + fileName;

    // Remove extension to create a clean base string.
    unique = stripExtension(fileName);

    rawAName = unique + "_raw";
    bgAName = unique + "_bgsub.tif";

    print("\\Clear");
    print("Processing: " + unique);
    print("Input file: " + fileName);
    print("Matched pattern: " + patternString);
    print("Rolling ball radius: " + rollingRadius + " px");

    // ---------- Open image ----------
    if (useBioFormats) {
        openImageNoDialog(pathA, rawAName);
    } else {
        open(pathA);
        rename(rawAName);
    }

    // ---------- Prepare background-subtracted image ----------
    selectWindow(rawAName);
    run("Duplicate...", "title=" + bgAName + " duplicate");
    selectWindow(bgAName);

    if (rollingRadius != 0) {
        run("Subtract Background...", "rolling=" + rollingRadius + " stack");
        bgApplied = "true";
    } else {
        bgApplied = "false";
    }

    // ---------- Optional save processed image ----------
    if (saveProcessedImages) {
        selectWindow(bgAName);
        saveAs("Tiff", outputDir + bgAName);
    }

    // ---------- Append log ----------
    File.append(
        csvEscape(unique) + "," +
        csvEscape(fileName) + "," +
        csvEscape(patternString) + "," +
        rollingRadius + "," +
        bgApplied + "," +
        csvEscape(bgAName) + "," +
        csvEscape("success") + "\n",
        logPath
    );

    // ---------- Close windows ----------
    closeIfOpen(rawAName);
    closeIfOpen(bgAName);

    run("Close All");
}

setBatchMode(false);
print("Done.");
print("Rolling-ball processing log saved to:");
print(logPath);


// ---------------- Helper functions ----------------

function closeIfOpen(title) {
    if (isOpen(title)) {
        selectWindow(title);
        close();
    }
}

function stripExtension(filename) {
    name = filename;

    // Remove longer extensions first.
    name = replace(name, ".ome.tiff", "");
    name = replace(name, ".ome.tif", "");
    name = replace(name, ".tiff", "");
    name = replace(name, ".tif", "");

    return name;
}

function csvEscape(s) {
    s = "" + s;
    s = replace(s, "\"", "\"\"");
    return "\"" + s + "\"";
}

function openImageNoDialog(path, newTitle) {
    // Non-interactive Bio-Formats import.
    // Good for OME-TIFF / Bio-Formats-readable files.
    run("Bio-Formats Importer",
        "open=[" + path + "] " +
        "color_mode=Default " +        
        "view=[Standard ImageJ] " +
        "stack_order=XYZCT"
    );

    rename(newTitle);
}