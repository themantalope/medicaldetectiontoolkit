import json
from metaflow import FlowSpec, IncludeFile, Parameter, step, JSONType
from filematcher import FileMatcher
from pprint import pprint


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

    base_dir = Parameter('base_dir')
    file_extension = Parameter('file_extension')
    out_meta_data_file = Parameter('out_meta_data_file')
    names_and_search_strings = Parameter('names_and_search_strings', type=JSONType, default='{}')
    main_name = Parameter('main_names')
    
    @step
    def start(self):
        fm = FileMatcher()
        # just get a generic list of unique ids
        names_dict = dict(self.names_and_search_strings)
        file_list = fm.collect_file_names(self.base_dir, [names_dict[self.main_name]])
        fm.find_unique_ids_in_file_list(file_list, name=names_dict[self.main_name])
        for k, v in names_dict.items():
            fl = fm.collect_file_names(self.base_dir, [names_dict[k]])
            fm.match_list_to_unique_ids(fl, v, k)
        pprint(fm.unique_dict)
        pass 

    @step
    def end(self):
        pass 


if __name__ == 'main':
    HCCPreProcessingFlow()