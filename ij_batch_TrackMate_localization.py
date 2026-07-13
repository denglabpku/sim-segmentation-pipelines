# ==========================================
# Batch TrackMate spot extraction from XML
# Zuhui Wang, 2026/03/07
#
# Purpose:
#   Load TrackMate settings from a saved XML template,
#   process images in batch, and export filtered spot
#   localizations to CSV.
#
# Inputs:
#   - template_xml_path : TrackMate XML template
#   - input_folder_path : folder containing images
#   - output_folder_path: folder for CSV output
#   - name_contain      : only process filenames containing this string
#   - exts_csv          : allowed file extensions
#   - n_threads         : number of CPU threads
#
# Output:
#   For each matching input image, save one CSV file with
#   the same basename as the image.
#
# Notes:
#	- Manually use TrackMate to provide the template XML first.
#   - Designed for single time-point z-stack images.
#   - Exports only spots remaining after TrackMate spot filtering.
#   - No tracking step is required.
# ==========================================

#@ String(label="Template TrackMate XML path", value="D:/path/to/template.xml") template_xml_path
#@ String(label="Input image folder", value="D:/path/to/images") input_folder_path
#@ String(label="Output folder", value="D:/path/to/output") output_folder_path
#@ String(label="Extensions (comma-separated)", value=".tif,.tiff,.nd2,.czi") exts_csv
#@ String(label="File name selection string", value="C1") name_contain
#@ Integer(label="Number of threads", value=4) n_threads
print('Hi')
import os
import sys

from java.io import File
from ij import IJ
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.io import TmXmlReader
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.io import CSVExporter

from java.io import BufferedWriter, FileWriter

def export_visible_spots_to_csv(csv_path, model):
    spots = model.getSpots()
    visible_spots = spots.iterable(True)   # True = visible spots only

    bw = BufferedWriter(FileWriter(csv_path))
    try:
        bw.write("ID,NAME,FRAME,X,Y,Z,RADIUS,QUALITY\n")
        for spot in visible_spots:
            sid = spot.ID()
            name = spot.getName()
            frame = spot.getFeature('FRAME')
            x = spot.getFeature('POSITION_X')
            y = spot.getFeature('POSITION_Y')
            z = spot.getFeature('POSITION_Z')
            radius = spot.getFeature('RADIUS')
            quality = spot.getFeature('QUALITY')

            bw.write("{},{},{},{},{},{},{},{}\n".format(
                sid, name, frame, x, y, z, radius, quality
            ))
    finally:
        bw.close()

reload(sys)
sys.setdefaultencoding('utf-8')

logger = Logger.IJ_LOGGER


def fail(msg):
    IJ.log("ERROR: " + msg)
    raise Exception(msg)


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def list_images(folder, exts):
    files = []
    for name in os.listdir(folder):
        full = os.path.join(folder, name)
        if not os.path.isfile(full):
            continue

        # Only process filenames containing "C1"
        if name_contain not in name:
            continue

        lower = name.lower()
        ok = False
        for ext in exts:
            if lower.endswith(ext.lower()):
                ok = True
                break
        if ok:
            files.append(full)

    files.sort()
    return files


def save_trackmate_xml(xml_path, model, settings, display_settings, log_text):
    out_file = File(xml_path)
    writer = TmXmlWriter(out_file, logger)
    writer.appendLog(log_text)
    writer.appendModel(model)
    writer.appendSettings(settings)
    if display_settings is not None:
        writer.appendDisplaySettings(display_settings)
        writer.appendGUIState('ConfigureViews')
    writer.writeToFile()


template_xml_file = File(template_xml_path)
if not template_xml_file.exists():
    fail("Template XML not found: " + template_xml_path)

reader = TmXmlReader(template_xml_file)
if not reader.isReadingOk():
    fail("Could not read template XML:\n" + reader.getErrorMessage())

template_display_settings = reader.getDisplaySettings()
template_log = reader.getLog()
if template_log is None:
    template_log = ""

ensure_dir(output_folder_path)

exts = [e.strip() for e in exts_csv.split(",") if e.strip()]
image_files = list_images(input_folder_path, exts)

if len(image_files) == 0:
    fail('No matching image files found containing "C1" in: ' + input_folder_path)

IJ.log("Found %d image(s) containing C1." % len(image_files))

n_ok = 0
n_fail = 0

for i, image_path in enumerate(image_files):
    image_name = os.path.basename(image_path)
    base_name = os.path.splitext(image_name)[0]

    IJ.log("")
    IJ.log("=== [%d/%d] Processing: %s ===" % (i + 1, len(image_files), image_name))

    imp = None
    try:
        imp = IJ.openImage(image_path)
        if imp is None:
            raise Exception("Could not open image: " + image_path)

        # Read settings from template XML and apply them to this image
        settings = reader.readSettings(imp)
        if settings is None:
            raise Exception("reader.readSettings(imp) returned None.")

        model = Model()
        model.setLogger(logger)

        tm = TrackMate(model, settings)
        tm.setNumThreads(n_threads)

        if not tm.checkInput():
            raise Exception("checkInput() failed:\n" + tm.getErrorMessage())

#        if not tm.process():
#            raise Exception("process() failed:\n" + tm.getErrorMessage())
#        if not tm.checkInput():
#			raise Exception("checkInput() failed:\n" + tm.getErrorMessage())

        if not tm.execDetection():
            raise Exception("Detection failed:\n" + tm.getErrorMessage())

        if not tm.computeSpotFeatures(True):
            raise Exception("Computing spot features failed:\n" + tm.getErrorMessage())

        if not tm.execInitialSpotFiltering():
            raise Exception("Initial spot filtering failed:\n" + tm.getErrorMessage())

        if not tm.execSpotFiltering(True):
            raise Exception("Spot filtering failed:\n" + tm.getErrorMessage())

        # Save TrackMate XML using same basename as input image
#        out_file_xml = os.path.join(output_folder_path, base_name + ".xml")
#
#        batch_log = (
#            "Batch-processed from template XML:\n"
#            + template_xml_path
#            + "\nInput image:\n"
#            + image_path
#            + "\n\nOriginal template log:\n"
#            + template_log
#        )

#        save_trackmate_xml(
#            out_file_xml,
#            model,
#            settings,
#            template_display_settings,
#            batch_log
#        )

        # Save spot localization CSV with same basename as input image
        out_file_csv = os.path.join(output_folder_path, base_name + ".csv")
#        only_visible = False
#        CSVExporter.exportSpots(out_file_csv, model, only_visible)
        export_visible_spots_to_csv(out_file_csv, model)

        IJ.log("OK: " + image_name)
#        IJ.log("Saved XML: " + out_file_xml)
        IJ.log("Saved CSV: " + out_file_csv)
        n_ok += 1

    except Exception as e:
        IJ.log("FAILED: " + image_name)
        IJ.log(str(e))
        n_fail += 1

    finally:
        if imp is not None:
            imp.close()

IJ.log("")
IJ.log("Batch finished.")
IJ.log("Succeeded: %d" % n_ok)
IJ.log("Failed: %d" % n_fail)
IJ.log("Output folder: " + output_folder_path)