import improc
import os 
import json
import SimpleITK as sitk
from metaflow import FlowSpec, Parameter, step, JSONType
from filematcher import FileMatcher
from pprint import pprint
from fuzzywuzzy import fuzz
from numpy import iterable, spacing, unique
from tqdm import tqdm



class HCCPreProcessingFlow(FlowSpec):
    
    """
    Steps in this flow:
    1) collect files
    2) organize files into a spreadsheet per patient
    3) if one name is the 'main' or central sequence, all other sequences should be aligned to it
    4) normalize data
    5) re-size data
    6) align data
    """
    
    base_dir = Parameter('base_dir', help='base dir to search')
    file_extension = Parameter('file_extension')
    out_meta_data_file = Parameter('out_meta_data_file')
    names_and_search_strings = Parameter('names_and_search_strings', type=JSONType, default='{}')
    main_name = Parameter('main_name')
    spacing = Parameter('spacing', type=JSONType, default='[1.0, 1.0, 1.0]')
    
    @step
    def start(self):
        fm = FileMatcher()
        # just get a generic list of unique ids
        names_dict = dict(self.names_and_search_strings)
        # pprint('names: {}'.format(names_dict))
        # pprint('names list: {}'.format([names_dict[self.main_name]]))
        file_list = fm.collect_file_names(self.base_dir, [names_dict[self.main_name]])
        fm.find_unique_ids_in_file_list(file_list, name=names_dict[self.main_name])
        for k, v in names_dict.items():
            fl = fm.collect_file_names(self.base_dir, [names_dict[k]])
            fm.match_list_to_unique_ids(fl, v, k)
        # pprint(fm.unique_dict)

        # now for a specific thing for this project
        for k in tqdm(fm.unique_dict.keys(),desc='matching arterial',total=len(fm.unique_dict.keys())):
            label_name = fm.unique_dict[k][self.main_name]
            # print(k)
            u_files = fm.collect_file_names(self.base_dir, [k], ext='.nii')
            u_files = list(set(u_files).difference({label_name}))
            # print(len(u_files))
            matches = sorted(
                [
                    (ufn, fuzz.partial_ratio(ufn.lower(), label_name.lower())) for ufn in u_files
                ], 
                key= lambda x: x[1]
            )
            best = matches[-1]
            # pprint(label_name)
            # pprint(best[0])
            fm.unique_dict[k]['arterial'] = best[0]
            # break

        print(fm.unique_dict)

        self.unique_dict = fm.unique_dict

        self.next(self.pre_resize)

    
    @step 
    def pre_resize(self):
        self.resize_list = []
        # need to resize everything
        for u, v in self.unique_dict.items():
            for n, f in v.items():
                self.resize_list.append((u, n, f))

        self.next(self.resize, foreach='resize_list') 

    @step
    def resize(self):
        unique, name, file = self.input 
        new_name = '{}{}'.format(name, '_processed')
        path, ext = os.path.splitext(file)
        new_outfile = '{}{}{}'.format(path, '_processed',ext)
        self.unique_dict[unique][new_name] = new_outfile
        self.new_file_updates = (unique, new_name, new_outfile)
        spacing = json.load(self.spacing)

        in_image = sitk.ReadImage(file)

        out_image = improc.resample(in_image, spacing)

        sitk.WriteImage(out_image, new_outfile)

        self.next(self.join_resize)
        

    @step
    def join_resize(self, inputs):

        for i in inputs:
            nfn = i.new_file_updates
            self.unique_dict[nfn[0]][nfn[1]] = nfn[2]

        self.merge_artifacts()
        self.next(self.pre_norm)
    
    
    @step
    def pre_norm(self):
        iterlist = []
        for k, v in self.unique_dict.items():
            iterlist.append((k, v['arterial_processed']))

        self.norm_list = iterlist
        self.next(self.norm, foreach='norm_list')

        
    @step
    def norm(self):
        unique, filename = self.input
        in_image = sitk.ReadImage(filename)
        out_image = improc.sitk_numpy_normalize(in_image)
        sitk.WriteImage(out_image, filename)
        self.new_file_updates = (unique, filename)
        self.next(self.norm_join)

    @step
    def norm_join(self):
        pass 


    @step 
    def pre_align(self):
        pass 

    @step
    def align(self):
        pass 

    @step
    def align_join(self):
        pass 



    
    @step
    def end(self):
        pass


if __name__ == '__main__':
    HCCPreProcessingFlow()
