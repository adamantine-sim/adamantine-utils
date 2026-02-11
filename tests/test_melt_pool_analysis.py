import os
import shutil
from adamantine_utils.melt_pool_analysis.melt_pool_analysis import melt_pool_statistics

def test_melt_pool_analysis():
    print("Testing melt-pool-analysis...")
    this_file_path = os.path.dirname(os.path.realpath(__file__))

    path_to_adamantine_files = this_file_path + "/test_vtu_files/"
    adamantine_filename = "solution"

    result = melt_pool_statistics(path_to_adamantine_files, adamantine_filename)

    print(result)

