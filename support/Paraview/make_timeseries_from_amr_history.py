#!/usr/bin/python

import sys
import os
import xml.etree.ElementTree as ET

def main(*args):
    roots = {}

    vtk_dir = args[0]
    for subdir_name in os.listdir(vtk_dir):
        if os.path.isdir(os.path.join(vtk_dir, subdir_name)):
            for file_name in os.listdir(os.path.join(vtk_dir, subdir_name)):
                if file_name.endswith('.pvtu'):
                    if file_name in roots:
                        root, collection = roots[file_name]
                    else:
                        root = ET.Element('VTKFile')
                        root.set("type", "Collection")
                        root.set("version", "0.1")
                        root.set("byte_order", "LittleEndian")
                        root.set("compressor", "vtkZLibDataCompressor")
                        collection = ET.SubElement(root, 'Collection')
                        roots[file_name] = (root, collection)
                    dataset = ET.SubElement(collection, 'DataSet')
                    dataset.set("timestep", subdir_name)
                    dataset.set("group", "")
                    dataset.set("part", "0")
                    dataset.set("file", os.path.join(subdir_name, file_name))
    for file_name, (root, _) in roots.iteritems():
        tree = ET.ElementTree(root)
        tree.write(os.path.join(vtk_dir, file_name.replace('.pvtu', '.pvd')), encoding='utf-8')

if __name__ == "__main__":
    main(*sys.argv[1:])
