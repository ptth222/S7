import numpy
import pandas
from multiprocessing import Pool
from functools import partial
import traceback

def get_loc(x, intensity_dataframe, lowest_index):
    ## Only search up through the next 500 m/z's. This drastically reduces the search space and execution time for large dataframes.
    low_index = x.name
    search_low = low_index - lowest_index
    try :
        return low_index + intensity_dataframe.iloc[search_low:search_low+500].index.get_loc(x["high_range"], method = "bfill")
    except KeyError:
        if search_low < (len(intensity_dataframe) - 499):
            print( "Max index in get_loc needs to be increased", low_index)
            raise Exception("Max index in get_loc needs to be increased", low_index)
        return (len(intensity_dataframe))
    
def mz_range_calcs(intensity_dataframe, low_index, high_index):
    try:
        num_of_signals_in_range = intensity_dataframe.iloc[low_index:high_index].count(axis=1).sum()
        if num_of_signals_in_range == 0:
            num_of_signals_in_range = 1
        
        mz_count = intensity_dataframe.iloc[low_index:high_index].count(axis=1)
        mz_count.index = range(len(mz_count))
        mean_intensity = intensity_dataframe.iloc[low_index:high_index].mean(axis=1).mean()/num_of_signals_in_range
        mean_mz = (pandas.Series(intensity_dataframe.index[low_index:high_index])*mz_count).sum()/num_of_signals_in_range
        weighted_mzs = numpy.repeat(intensity_dataframe.index[low_index:high_index], mz_count)
        
        return (num_of_signals_in_range, 
                intensity_dataframe.iloc[low_index:high_index].sum(axis=1).sum()/num_of_signals_in_range,
                numpy.nanstd((intensity_dataframe.iloc[low_index:high_index].values/num_of_signals_in_range), ddof=1)/mean_intensity*100,
                numpy.median(weighted_mzs), 
                numpy.nanstd(weighted_mzs, ddof=1)/mean_mz*100)
    except Exception as e:
        print( "Caught exception in worker thread mz_range_calcs:", num_of_signals_in_range)
        traceback.print_exc()
        print()
        raise e


def parallelize_dataframe(df, func, num_of_cores, intensity_dataframe):
    chunk_size = len(df) // num_of_cores + 1
    indices = range(chunk_size, (len(df) // chunk_size)*chunk_size, chunk_size)
    if len(indices) == 1:
        indices = 1
    df_split = numpy.split(df, indices)
    
    intensity_df_split = []
    for frame in df_split:
        intensity_df_split.append(intensity_dataframe.iloc[frame.iloc[0].name:(frame.iloc[-1].name+500)])
    pool = Pool(num_of_cores)
    df = pandas.concat(pool.map(func, zip(df_split, intensity_df_split)))
    pool.close()
    pool.join()
    return df

def parallel_range_calcs(df_intensity_df_tuple):
    try:
        data = df_intensity_df_tuple[0]
        intensity_dataframe = df_intensity_df_tuple[1]
    
        temp = \
        pandas.DataFrame(
                list(
                        data.apply(
                                lambda x: 
                                    mz_range_calcs(
                                            intensity_dataframe,
                                            int(x["low_index"])-data.index[0], 
                                            int(x["high_index"])-data.index[0]
                                            ), 
                                    axis = 1
                                    )
                                    ), 
                            
                columns=["Number of Signals in m/z Range", 
                         "(Sum of Intensities)/(Number of Scans)", 
                         "RSTD of Intensities/(Number of Scans)", 
                         "Median m/z in m/z Range", 
                         "m/z RSTD in m/z Range"]
                )
        return temp
    except Exception as e:
        print( "Caught exception in worker thread parallel_range_calcs:")
        traceback.print_exc()
        print()
        raise e





def parallel_get_loc(df_intensity_df_tuple):
    data = df_intensity_df_tuple[0]
    intensity_dataframe = df_intensity_df_tuple[1]
    data["high_index"] = data.apply(lambda x: get_loc(intensity_dataframe = intensity_dataframe, x = x, lowest_index = data.index[0]), axis=1)  
    return data


def save_writer(writer):
    writer.save()

def parallel_save(writers, num_of_cores):
    pool = Pool(num_of_cores)
    pool.map(save_writer, writers)



def f(x):
    return x*x

def test():
    p = Pool(2)
    print( p.map(f, [1, 2, 3, 4, 5, 6]))
    
def df_square (df):
    return df*2
