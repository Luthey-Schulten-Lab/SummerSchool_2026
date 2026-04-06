"""
Conveinent mathmatical functions useful in WCM analysis
"""
import numpy as np

def cal_CV(array):
    """
    Calculate the coefficient of variation of array
    """
    
    CV = np.std(array,axis=None, ddof=1)/np.mean(array)

    return CV

def log_0(arr):
    """
    skip negative when taking log
    """

    return np.where(arr > 0, np.log(arr), np.nan)


def divide_replace(numerator, denominator, value):
    """
    Divide numerator by denominator and set nan and inf to value
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        result = numerator / denominator
        # Replace NaN and Inf values with finite value
        result[np.isnan(result) | np.isinf(result)] = value

    return result

def get_min_median_mean_max(list):
    """
    return a list of mean, median, mean and max of a 1D numpy array
    """
    return np.array([np.min(list), np.median(list), np.mean(list), np.max(list)])

def get_doubling_moments_value(doubling_times, array):
    """
    Input: 
    doubling_times, a list,  doubling times for each replicate
    array, nD (n=2 or 3)numpy array where the last dimesion is replicates
    
    Return: 
    value, (n-1)D numpy array by extracting the value at time moment equals to doubling time 
    """

    doubling_times = [int(time) for time in doubling_times]
    
    N_reps = np.shape(array)[-1]
    
    # print('N_reps, {0}'.format(N_reps))

    if len(doubling_times) != N_reps:
        print("WARNING: Inputs Dimensions Don't Match")
        return None

    if array.ndim == 2:
        value = array[doubling_times,np.arange(N_reps)]
    elif array.ndim == 3:
        value = array[:,doubling_times, np.arange(N_reps)]
    else:
        print("WARNING: Dimension of Input is {0}, Not valid array for Function get_doubling_moments_value".format(array.ndim))
        return None

    return value

def get_mean_upto_moments(doubling_times, array, start=0):
    """
    Input: 
    doubling_times, a list,  doubling times for each replicate
    array, nD (n=2 or 3)numpy array where the last dimesion is replicates
    start: starting time point, in seconds

    Return: 
    value, (n-1)D numpy array by calculating the time average up to the doubling times
    """


    doubling_times = [int(time) for time in doubling_times]
    
    N_reps = np.shape(array)[-1]
    
    # print('N_reps, {0}'.format(N_reps))
    def cal(twoDarray, doubling_times):

        truncated = [twoDarray[start:doubling_times[i],i] for i in range(len(doubling_times))]

        means = np.array([np.mean(_) for _ in truncated])

        return means

    if len(doubling_times) != N_reps:
        print("WARNING: Inputs Dimensions Don't Match")
        return None

    if array.ndim == 2: # time, replicates
        value = cal(array, doubling_times)

    elif array.ndim == 3: # species, time, replicates
        value = np.zeros((array.shape[0], array.shape[1]))
        for i in range(array.shape[0]):
            value[i,:] = cal(array[i,:,:], doubling_times)

    else:
        print("WARNING: Dimension of Input is {0}, Not valid array for Function get_doubling_moments_value".format(array.ndim))
        return None

    return value

def replace_symbols(s, replacements):
    """
    
    Convert strings to raw string while converting the letters
    """
    import re
    # replacements = {
    #     'a': r'\alpha',
    #     'b': r'\beta',
    #     'g': r'\gamma',
    #     'd': r'\delta',
    #     'e': r'\epsilon',
    #     'A': ' a',
    #     'B': ' b',
    #     'C': ' c'  # Convert uppercase 'C' to lowercase
    # }

    # Dynamically create regex pattern from the dictionary keys
    pattern = r'[' + ''.join(re.escape(key) for key in replacements.keys()) + r']'

    # Use regex to replace each character with its corresponding LaTeX string
    def replace_match(match):
        char = match.group(0)
        return replacements.get(char, char)  # Replace if in dict, else keep as is

    # Apply the replacements
    transformed_s = re.sub(pattern, replace_match, s)
    
    # Add 'r' in front for LaTeX rendering
    return r'$' + transformed_s + '$'

