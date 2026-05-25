import xml.etree.ElementTree as ET


class XMLDataExtractor:
    def __init__(self, xml_file_path):
        self.xml_file_path = xml_file_path
        self.tree = ET.parse(self.xml_file_path)
        self.root = self.tree.getroot()
    

    def get_box_size(self):
        box_size: list[float] = [float(self.root.find('.//box').attrib[i]) for i in ['lx', 'ly', 'lz']]
        print(f"Your box size: a = {box_size[0]}, b = {box_size[1]}, c = {box_size[2]}.")
        return box_size
