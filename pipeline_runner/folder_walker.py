import json
import os

class FolderWalker:
    """A walker that walk through the folder tree, and convert it to JSON data which is required by jsTree.js.

    Attributes
    ----------
    folder_path : str
        Parent directory of the folder to be parsed.
    folder : str
        Name of the folder to be parsed.
    data : list
        A list of dict represents the folder structure.
    """
    def __init__(self, folder):
        self.folder_path = os.path.dirname(os.path.abspath(folder))
        self.folder = os.path.basename(folder)
        self.data = self.get_folder()

    def get_folder(self):
        """Walk through the whole folder structure, construct a list of dict represents the folder structure.

        Returns
        -------
        data : list
            A list of dict represents the folder structure.
        """
        current_dir = os.getcwd()
        os.chdir(self.folder_path)
        
        data = []
        data.append({'parent' : '#',
                     'id' : self.folder,
                     'text' : self.folder})
        
        for root, dirs, files in os.walk(self.folder):
            for d in dirs:
                x = {
                    'parent' : root,
                    'id' : os.path.join(root, d),
                    'text' : d,
                }    
                data.append(x)
            for f in files:
                x = {
                    'parent' : root,
                    'id' : os.path.join(root, f),
                    'text' : f,
                    'icon' : 'far fa-file',
                    'a_attr' : { 'href' : os.path.join(root, f)}
                }
                data.append(x)

        os.chdir(current_dir)
        return data

    def to_json(self, fname):
        """Write the data to JSON file.

        Returns
        -------
        None
        """
        with open(fname, 'w') as f:
            json.dump(self.data, f)
