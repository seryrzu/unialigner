import pandas as pd
import numpy as np


def get_bridges(bridges_csv_path):
    bridges_df = pd.read_csv(bridges_csv_path, header=1)
    bridges_df["na"] = [np.nan] * len(bridges_df)
    print(bridges_df.columns)
    bridges_X = np.ravel(bridges_df[["start_1", "end_1", "na"]])
    bridges_Y = np.ravel(bridges_df[["start_2", "end_2", "na"]])
    return bridges_X, bridges_Y
