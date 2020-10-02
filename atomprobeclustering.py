import numpy as np
import os
import pandas
import plotly.graph_objects as go
import plotly.express as px
from plotly.colors import n_colors
import plotly.figure_factory as ff
from scipy import stats
from typing import List, Optional, Dict

# constant values and lists useful for data processing
constant_cluster_stats_columns = ["X", "Y", "Z", "Unranged", "r_gyration"]

# TODO add all atomic densities for every element - is there any external library with the densities?
atomic_densities = {
    "Mo" : 0.01558,
    "Fe" : 0.01177,
    "Ni" : 0.01094,
    "Cu" : 0.01181,
    "Cr" : 0.01201,
    "Si" : 0.02003,
    "P" : 0.02826,
    "O" : 0.02883,
    "C" : 0.00878,
    "Mn" : 0.01221,
    "Ga" : 0.01960,
    "MoN" : 0.03806,
    "C2" :0.01757, 
    "C3" :0.02635 
}

# TODO is there a way to import all IVAS colours?
core_ion_colors = {
    "Cu": "orange",
    "Ni": "green",
    "Mn": "yellow",
    "P": "pink",
    "Si": "grey"
}

# colors = n_colors('rgb(5, 200, 200)', 'rgb(200, 10, 10)', 6, colortype='rgb')

def chart_colors(n:int) -> list:
    """Return iterable list of n colors

    Args:
        n (int): how many colors do you need?

    Returns:
        list: n rgb colors
    """
    colors = n_colors('rgb(5, 200, 200)', 'rgb(200, 10, 10)', n, colortype='rgb')
    return colors


string_vector = List[str]
float_vector = List[float]
integer_vector = List[int]
dict_with_pandas_dfs = Dict[int, pandas.DataFrame]
# read and process the data and combine into one pandas dataframe

def read_cluster_file(
    cluster_stats_file:str,
    core_ions:list,
    exclude_ions:list=[],
    sample_title:str="",
    output_folder:str="",
    n_min:int=1
) -> pandas.DataFrame:
    
    if len(sample_title) != 0:
        sample_title = cluster_stats_file.replace("_cluster-stats.txt", "").replace(".txt", "")
    
    # open the cluster stats file and create pandas dataframe
    local_df = pandas.read_csv(f"{output_folder}{cluster_stats_file}", sep="\t")
    # TODO: change to full file path?

    # insert all the corrections in here
    # exclude (matrix) ions by deleting the column
    local_df = local_df.drop(exclude_ions, axis=1)
    
    # # decompose columns if needed (this block would only give it correct number of ions needed for total cluster size)
    # for col in local_df.columns:
    #     if col not in constant_cluster_stats_columns:
    #         local_decomposed_ions = 1
    #         for letter in col:
    #             if letter.isupper():    # every element starts with capital letter
    #                 local_decomposed_ions += 1
    #             elif unicode(letter).isnumeric(): # some elements are multiplied
    #                 local_decomposed_ions += (int(letter) - 1)
    #         # multiply the number of these ions so that the number corresponds to real number or atoms in each cluster
    #         local_df[col] = local_df[col] * local_decomposed_ions

    # add total number of core ions and create a column with cluster size    
    # local_df["Cluster Total Size"] = local_df.iloc[:, 3:-2].sum(axis=1)
    local_df["Cluster Size"] = local_df[core_ions].sum(axis=1)
    
    # filter by n_min
    local_df = local_df[local_df["Cluster Size"] >= n_min]
    
    # sort it by size
    local_df = local_df.sort_values(by="Cluster Size")
    
    return local_df

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def decompose_columns(
    local_df:pandas.DataFrame, 
    skip_columns:string_vector=constant_cluster_stats_columns
) -> pandas.DataFrame:    
    """ work in progress, don't rely on it

    Args:
        local_df (pandas.DataFrame): [description]
        skip_columns (string_vector, optional): [description]. Defaults to constant_cluster_stats_columns.

    Returns:
        pandas.DataFrame: [description]
    """
    # make another decomposition function which assigns ions from the block above to the correct columns and empties (unde)composed columns
    for col in local_df.columns:
        if col not in skip_columns:
            local_elements = {}
            local_element = col[0]  # first letter from the column
            for letter in col[1:]:
                if letter.isnumeric():
                    local_elements[local_element] = int(letter)
                    local_element = ""                    
                elif letter.isupper():
                    local_elements[local_element] = 1
                    local_element = letter
                elif letter.islower():
                    local_element += letter
                else:
                    print("Something wrong with decomposition.")
            
            # assign the last element
            local_elements[local_element] = 1 
            local_element = ""

            # assign these ions across other columns, create new if needed
            for key, value in local_elements.items():
                if key not in local_df.columns:
                    local_df[key] = local_df[col] * value
                else:
                    local_df[key] += local_df[col] * value
    return local_df

def cluster_summary(
    cluster_stats:pandas.DataFrame,
    unclustered_stats:pandas.DataFrame,
    sample_title:str,
    # n_min:int=1,
    eta:float=0.52, # assuming it's LEAP 5000 XR with 52% detection efficiency
    matrix_ions:string_vector=["Fe"]
) -> pandas.DataFrame:
    
    """Function that reads file path to cluster-stats file from posgen cluster search
    and returns volume fraction, number density, average radius, all with or without matrix ions,
    and their errors. Uses cluster-stats.txt and unclustered.txt files.

    Raises:
        FileNotFoundError: [description]
        FileNotFoundError: [description]
        FileNotFoundError: [description]

    """

    # TODO change with_Fe into with_matrix, it's supposed to be used with non-ferritic materials too
    # ask about the matrix element to create column names

    volume_fractions_with_matrix_atoms = []
    volume_fractions_no_matrix_atoms = []
    number_densities = []
    average_radii = []
    average_radii_error = []
    volume_fractions_with_matrix_atoms_simple = []
    volume_fractions_no_matrix_atoms_simple = []
    volume_fraction_error = []
    volumes = []

    # start iterating here

    # find names of the files
    # cluster_file = str(cluster_file_path.split(sep="/")[-1])
    # output_folder = cluster_file_path.replace(cluster_file, "") # find the directory where the cluster file is
    # unclustered_file = cluster_file.replace("_cluster-stats.txt", "_unclustered-stats.txt")

    # open both files in pandas
    # cluster_stats = pandas.read_csv(output_folder+cluster_file, sep="\t")
    # unclustered_stats = pandas.read_csv(output_folder+unclustered_file, sep="\t")

    # calculate volumes
    # adjust Fm column - bad overlap solver for peak 29 Da
    # cluster_df["Fe"] = cluster_df["Fe"] + (cluster_df["Fm"] * 0.1)
    # cluster_df["Ni"] = cluster_df["Ni"] + (cluster_df["Fm"] * 0.9)

    # filter clusters below N_min threshold
    # n_min = n_mins[i] # TODO delete
    cluster_df_processed = cluster_stats
    # cluster_df_processed["Cluster Size"] = cluster_df_processed.iloc[:, 3:-2].sum(axis=1)
    # cluster_df_processed = cluster_df_processed[cluster_df_processed["Cluster Size"] >= n_min]

    # clusters with Fe
    cluster_volume_with_matrix_atoms = 0
    cluster_volume_no_matrix_atoms = 0
    unclustered_volume_with_matrix_atoms = 0

    # longer method
    for ion in atomic_densities:
        try:
            cluster_volume_with_matrix_atoms += cluster_stats[ion].sum() * atomic_densities[ion] # units: atom * nm3/atom = nm3
            if ion in matrix_ions:
                cluster_volume_no_matrix_atoms += cluster_df_processed[ion].sum() * atomic_densities[ion] # units: atom * nm3/atom = nm3, don't include Fe
            unclustered_volume_with_matrix_atoms += unclustered_stats[ion].sum() * atomic_densities[ion] # units: atom * nm3/atom = nm3
        except KeyError:
            print("{} ion column not found in {}".format(ion, cluster_stats))
    tip_volume = cluster_volume_with_matrix_atoms + unclustered_volume_with_matrix_atoms
    this_volume_fraction_no_matrix_atoms = cluster_volume_no_matrix_atoms / tip_volume
    this_volume_fraction_with_matrix_atoms = cluster_volume_with_matrix_atoms / tip_volume
    
    # simplified volume fraction cluster ion count / total ion count
    cluster_count_with_matrix_atoms = 0
    cluster_count_no_matrix_atoms = 0
    unclustered_count_with_matrix_atoms = 0

    for ion in atomic_densities:
        try:
            cluster_count_with_matrix_atoms += cluster_stats[ion].sum()
            if ion in matrix_ions:
                cluster_count_no_matrix_atoms += cluster_df_processed[ion].sum()
            unclustered_count_with_matrix_atoms += unclustered_stats[ion].sum()
        except KeyError:
            print("{} ion column not found in {}".format(ion, cluster_stats))
    tip_ion_count = cluster_count_with_matrix_atoms + unclustered_count_with_matrix_atoms
    this_volume_fraction_no_matrix_atoms_simple = cluster_count_no_matrix_atoms / tip_ion_count
    this_volume_fraction_with_matrix_atoms_simple = cluster_count_with_matrix_atoms / tip_ion_count    

    # estimating error of volume fraction
    Nt = float(tip_ion_count)
    Nc = float(this_volume_fraction_no_matrix_atoms_simple)
    if (Nt == 0) or (Nc == 0):
        this_volume_fraction_error = (Nt**(-2)*(Nc**(-1/2))+((Nc/Nt)**2 * (1/Nt)))**0.5
    else:
        this_volume_fraction_error = None

    # number_density = cluster count / (LEAP efficiency * atomic volume of the tip)
    this_number_density = len(cluster_stats.index) / (eta * tip_volume)

    # average radius and st dev
    this_average_radius = cluster_df_processed["r_gyration"].mean()
    this_average_radius_stdev = cluster_df_processed["r_gyration"].std()

    # export to the lists
    volume_fractions_no_matrix_atoms.append(this_volume_fraction_no_matrix_atoms)
    volume_fractions_with_matrix_atoms.append(this_volume_fraction_with_matrix_atoms)
    number_densities.append(this_number_density)
    average_radii.append(this_average_radius)
    average_radii_error.append(this_average_radius_stdev)
    volume_fractions_no_matrix_atoms_simple.append(this_volume_fraction_no_matrix_atoms_simple)
    volume_fractions_with_matrix_atoms_simple.append(this_volume_fraction_with_matrix_atoms_simple)
    volume_fraction_error.append(this_volume_fraction_error)
    fixed_volume = eta * tip_volume
    volumes.append(fixed_volume)

    # TODO: add error values based on sigma = sqrt(N) with correct propagation

    # print("\n")
    # for sample_title, vf, nd, r in zip(sample_titles, volume_fractions_no_matrix_atoms, number_densities, average_radii):
    #     print("{}:\t {:.3f}% | {:.3f} nm3 | {:.3f} nm".format(sample_title, vf*100, nd, r))


    # create pandas dataframe and populate the columns
    df_comparison = pandas.DataFrame()
    # df_comparison["Ageing Time"] = ages
    df_comparison["Volume Fraction"] = volume_fractions_no_matrix_atoms
    df_comparison["Number Density"] = number_densities
    df_comparison["Average Radius"] = average_radii
    # df_comparison["Sample"] = [file_name.replace("_o7_cluster-stats.txt", "") for file_name in cluster_stats_files]
    df_comparison["Volume Fraction (count)"] = volume_fractions_no_matrix_atoms_simple
    df_comparison["Volume Fraction Error"] = volume_fraction_error
    df_comparison["Volume [nm3]"] = volumes
    df_comparison["Sample"] = sample_title

    return df_comparison

def combine_cluster_summaries(*summary_dataframes:pandas.DataFrame) -> pandas.DataFrame:
    """combines multiple pandas dataframes from cluster summaries

    Returns:
        pandas.DataFrame: cluster summary for multiple samples
    """

    dataframes = []
    for df in summary_dataframes:
        dataframes.append(df)
    
    return pandas.concat(dataframes)

def plot_column_chart_cluster_composition(
    cluster_dataframe:pandas.DataFrame, 
    core_ions:string_vector, 
    bulk_ions:string_vector,
    plot_bulk_ions:bool = True,
    n_min:int = 1,
    eta:float = 1.0,
    sample_title:str = ""
) -> None:
    """Plot stacked column chart with composition of core ions in each cluster.
    Clusters are sorted by size and corrected for detection efficiency (eta).

    Args:
        cluster_dataframe (pandas.DataFrame): [description]
        core_ions (string_vector): [description]
        bulk_ions (string_vector): [description]
        plot_bulk_ions (bool, optional): [description]. Defaults to True.
        n_min (int, optional): [description]. Defaults to 1.
        eta (float, optional): [description]. Defaults to 1.0.
        sample_title (str, optional): [description]. Defaults to "".
    """
    
    # for i, cluster_stats_file in enumerate(cluster_stats_files):
    # copy dataframe
    local_df = cluster_dataframe
    
    # create a new figure
    fig = go.Figure()
    
    # # open the cluster stats file and create pandas dataframe
    # sample_title = cluster_stats_file.replace("_cluster-stats.txt", "")
    # local_df = pandas.read_csv(f"{output_folder}{cluster_stats_file}", sep="\t")
    
    # # insert all the corrections in here
    # exclude (matrix) ions by deleting the column
    local_df = local_df.drop(bulk_ions, axis=1)
    # # lower the amount of ions for each element that has overlap with matrix ions
    # for overlap_ion in overlap_ions:
    #     local_df[overlap_ion] = round(local_df[overlap_ion] * overlap_ions[overlap_ion])
    # # transfer the overlapping ions into "original" elements
    # local_df["Ni"] = local_df["Ni"] + local_df["Fm"]
    # # delete the overlaping ions by multiplying by 0 to keep the total amount of ions constant
    # local_df["Fm"] = local_df["Fm"] * 0
    
    # # add total ions in cluster value
    # local_df["Cluster Size"] = local_df.iloc[:, 3:-2].sum(axis=1)
    
    # change n_min
    local_df = local_df[local_df["Cluster Size"] >= n_min]
    
    # sort it by size
    local_df = local_df.sort_values(by="Cluster Size")
    
    previous_y_values = []
    # create x-axis based on cluster size 
    x_values = np.arange(len(local_df["Cluster Size"]))
    # consider detection efficiency and divide cluster sizes by this number
    x_values = x_values / eta

    # iterate through rows and find their core ions and total ions values
    for i2, core_ion in enumerate(core_ions):
        y_values = local_df[core_ion] / local_df["Cluster Size"] * 100
        
        # # for each cluster size find all clusters and calculate mean core ion content with stdev
        # for i3, x in enumerate(iter(x_values)):
        #     all_for_one_size = local_df[local_df["Cluster Size"] == x][core_ion]
        #     y_values[i3] = all_for_one_size.mean() / x * 100
        #     y_std[i3] = all_for_one_size.std() / x * 100
        #     # TODO: think about how to plot st dev later on

        fig.add_trace(go.Bar(
            x=x_values, 
            y=y_values, 
            name=core_ion,
            marker_color=core_ion_colors[core_ion]
        ))

    # # assign for later use    
    # dfs[i] = local_df

    # reopen pandas dataframe to get the Fe content and plot a line
    if plot_bulk_ions:
        # local_df_Fe = pandas.read_csv(f"{output_folder}{cluster_stats_file}", sep="\t")
        df_with_bulk_ions = cluster_dataframe
                                
        # # transfer the overlapping ions into "original" elements
        # local_df_Fe["Fe"] = local_df_Fe["Fe"] + (local_df_Fe["Fm"] * 0.1)
        # # delete the overlaping ions by multiplying by 0 to keep the total amount of ions constant
        # local_df_Fe["Ni"] = local_df_Fe["Ni"] + (local_df_Fe["Fm"] * 0.9)

        # add total ions in cluster value
        df_with_bulk_ions["Total Ions"] = df_with_bulk_ions.iloc[:, 3:-2].sum(axis=1)

        # change n_min
        df_with_bulk_ions = df_with_bulk_ions[df_with_bulk_ions["Total Ions"] >= n_min]

        # sort it by size
        df_with_bulk_ions = df_with_bulk_ions.sort_values(by="Cluster Size")

        # create column with sum of all bulk ions and find the percentage to total ions in each cluster
        df_with_bulk_ions["Total Bulk Ions"] = df_with_bulk_ions[bulk_ions].sum(axis=1)
        bulk_ions_concentration = df_with_bulk_ions["Total Bulk Ions"] / df_with_bulk_ions["Total Ions"] * 100

        # TODO: data smoothing?
        # plt.plot(x_values, Fe_concentration, 'r-', label="Fe content", alpha=0.7)
        bulk_line_name = ""
        for ion in bulk_ions:
            bulk_line_name += ion + "+"
        bulk_line_name = bulk_line_name[:-1]

        fig.add_trace(go.Line(x=x_values, y=bulk_ions_concentration, name=bulk_line_name, opacity=0.7, line_color="red"))
        
    # adjust the plot
    fig.update_layout(
        barmode="stack", 
        title_text=sample_title
    )

    # Set custom x-axis labels
    how_many_ticks = 15
    print(len(x_values))
    interval = round(len(x_values) / how_many_ticks)
    
    ticktext_data = local_df["Cluster Size"][::interval]
    ticktvals_data = x_values[::interval]
    
    fig.update_xaxes(
        ticktext=ticktext_data,
        tickvals=ticktvals_data,
        title_text="Cluster Size in Ions"
    )

    fig.update_yaxes(
        title_text="Concentration in at. %"
    )

    fig.show()

    #     # find the sum of ternary ions in each cluster    
    # local_df["Relative Cluster Size"] = local_df["Ni"] + local_df["Cu"] + local_df["Mn"] + local_df["Si"]
    
    # # once you have sum of ternary ions in each cluster, find relative concentration of each core ion
    # local_df["Ni+Mn"] = local_df["Ni"] + local_df["Mn"]
    # # local_df["Relative Si"] = local_df["Si"] / relative_cluster_size
    # # local_df["Relative Cu"] = local_df["Cu"] / relative_cluster_size

    # # create column with negative z location 
    # local_df["Z Axis"] = local_df["Z"] * -1

    return None

def plot_all_clusters_ternary_diagram(
    cluster_dataframe:pandas.DataFrame,
    ternary_elements:string_vector,
    n_min:int = 1,
    sample_title:str = ""
) -> None:
    """plots ternary diagram with composition of all clusters

    Args:
        cluster_dataframe (pandas.DataFrame): [description]
        ternary_elements (string_vector): [description]
        n_min (int, optional): [description]. Defaults to 1.
        sample_title (str, optional): [description]. Defaults to "".

    Returns:
        [type]: [description]
    """

    # for i, cluster_stats_file in enumerate(cluster_stats_files):
    # create a new figure
    # fig = go.Figure()
    # n_min = n_mins[i]
    # open the cluster stats file and create pandas dataframe
    # sample_title = cluster_stats_file.replace("_cluster-stats.txt", "")
    local_df = cluster_dataframe
     
    # change n_min
    local_df = local_df[local_df["Cluster Size"] >= n_min]
    
    # sort it by size
    # local_df = local_df.sort_values(by="Cluster Size")
    
    # previous_y_values = []
    # # iterate through rows and find their core ions and total ions values
    # x_values = np.arange(len(local_df["Cluster Size"]))
    # for i2, core_ion in enumerate(core_ions):
    #     y_values = local_df[core_ion] / local_df["Cluster Size"] * 100
        
    #     # fig.add_trace(go.Bar(
    #     #     x=x_values, 
    #     #     y=y_values, 
    #     #     name=core_ion,
    #     #     marker_color=core_ion_colors[core_ion]
    #     # ))
    
    # once you have sum of ternary ions in each cluster, find relative concentration of each core ion
    # local_df["Mn+Si"] = local_df["Si"] + local_df["Mn"]

    # find ternary elements, and create sums if required
    for ternary_element in ternary_elements:
        if "+" in ternary_element:
            # create a list of elements and sum them up
            sub_ternary_elements = ternary_element.split(sep="+")
            local_df[ternary_element] = local_df[sub_ternary_elements].sum(axis=1)


    # find the sum of ternary ions in each cluster    
    # local_df["Number of ternary ions"] = local_df[ternary_elements].sum(axis=1)

    # create column with negative z location 
    # local_df["Z Axis"] = local_df["Z"] * -1

    # plotting using plotly ternary diagrams
    
    fig2 = px.scatter_ternary(
        local_df, a=local_df[0], b=local_df[1], c=local_df[2], 
        hover_name="Cluster Size", # change later to show cluster ID
        size="Cluster Size", size_max=10,
        color="Cluster Size", color_continuous_scale=px.colors.sequential.Viridis,
        opacity=0.6
        # marker=dict(
        #     color='LightSkyBlue',
        #     size='Relative Cluster Size',
        #     size_max = 15,
        #     opacity=0.5,
        #     line=dict(
        #         color='MediumPurple',
        #         width=5
        #     )
        # )
    )

    fig2.update_layout(title_text=sample_title)
    fig2.show()

    return None

def plot_comparison_binned_compositions(
    cluster_dataframes:dict_with_pandas_dfs,
    y_axis:str,
    x_axis:str,
    bins:float_vector,
    sample_titles:Optional[string_vector] = None
) -> None:

    # copy dfs, need to be a list containing pandas.DataFrames 
    dfs = cluster_dataframes

    # TODO check if plotting function don't have to return fig to easily save them
    # normalise x-axis so that cluster size considers eta for each test
    fig = go.Figure()
    # y_axis = "Cu"
    # bins = np.arange(0, 500, 50)

    # set colors for the samples
    bar_colors = chart_colors(len(dfs))

    for i, df in dfs:
        # multiple columns can be added (i.e. concentration of 2+ elements in clusters)
        if "+" in y_axis:
            # create a list of elements and sum them up
            sub_columns = y_axis.split(sep="+")
            dfs[y_axis] = dfs[sub_columns].sum(axis=1)

        # find bin stats to smooth the graphs and generalise the data
        local_x_values = df[x_axis]
        local_comp = df[y_axis] * 100 / df["Cluster Size"] # change into percentage
        # local_x_values = local_cs / etas[i] # cluster size / LEAP efficiency

        # iterate through every bin and find mean and std of ion
        y_axis_values = np.zeros(len(bins)-1)
        y_axis_error = np.zeros(len(bins)-1)

        for i2, bin in enumerate(list(bins[1:])):
            not_lower_than = local_x_values > bins[i2]
            not_higher_than = local_x_values < bins[i2+1]
            mask = not_higher_than & not_lower_than
            bin_ion_mean = np.mean(local_comp[mask])
            bin_ion_std = np.std(local_comp[mask])
            try:
                bin_error = bin_ion_std / (len(local_comp[mask])**0.5) # SEM
            except ZeroDivisionError:
                bin_error = 0
            y_axis_values[i2] = bin_ion_mean
            y_axis_error[i2] = bin_error
        
        fig.add_trace(go.Bar(
            x=bins,
            y=y_axis_values, 
            name=sample_titles[i], 
            marker_color=bar_colors[i],
            # mode="markers",
            # opacity=1,
            # marker_line_width=20,
            error_y=dict(
                type='data', # value of error bar given in data coordinates
                # array=Cu_comp_std/(Cu_count**0.5), # error = SEM based on sigma / sqrt(N)
                array=y_axis_error, # error = sigma
                visible=True,
                width=1
                )
        ))
        # break

    # Overlay both histograms
    fig.update_layout(
        yaxis_title=f"{y_axis} mean composition [at.%]",
        title_text=x_axis,
        barmode='group' # group together boxes of the different traces for each value of x
    )

    # set up the x axis
    ticktext = []
    tickvals = []
    for i in range(len(bins[1:])):
        text = "{}-{}".format(bins[i], bins[i+1])
        ticktext.append(text)
        tickvals.append(bins[i]) 

    fig.update_xaxes(
        ticktext=ticktext,
        tickvals=tickvals,
        title_text=x_axis
    )

    # TODO check if that works in non-jupyter files and if you can save it later
    # TODO i think it'll require ' return fig ' instead 
    fig.show()

    return None

# TODO create function to combine all cluster dataframes into one dictionary

def plot_comparison_cluster_size_distributions(
    cluster_dataframes:dict_with_pandas_dfs,
    sample_titles:string_vector,
    x_axis:str="r_gyration"
) -> None:

    # Group data together
    hist_data = []

    group_labels = []

    for i, sample in enumerate(sample_titles):
        group_labels.append(sample)
        hist_data.append(cluster_dataframes[i][x_axis])

    # Create distplot with custom bin_size
    fig = ff.create_distplot(
        hist_data, 
        group_labels, 
        # bin_size=0.01
        show_hist=False
    )
    fig.update_layout(
        xaxis_title=x_axis,
        title_text="Cluster Size Distribution"
        # barmode='group' # group together boxes of the different traces for each value of x
    )
    fig.show()

    return None

def plot_comparison_binned_number_densities(
    cluster_dataframes:dict_with_pandas_dfs,
    bins:float_vector,
    cluster_summaries:pandas.DataFrame,
    x_axis:str="r_gyration",
    sample_titles:Optional[string_vector] = None
) -> None:

    # plot spatial distribution with number density of all clusters on y-axis and r_gyration on x-axis
    fig = go.Figure()

    bins = np.arange(0.2, 1.8, 0.2)
    df_comparison = cluster_summaries
    dfs = cluster_dataframes

    # set up colors for samples
    bar_colors = chart_colors(len(cluster_dataframes))

    for i, sample in enumerate(sample_titles):
        
        # find bin stats to smooth the graphs and generalise the data
        # local_cs = dfs[i]["Cluster Size"]
        # local_comp = dfs[i][ion] / dfs[i]["Cluster Size"]
        # local_cs_normalised = local_cs / etas[i] # cluster size / LEAP efficiency
        local_x_values = dfs[i][x_axis]

        # iterate through every bin and find mean and std of ion
        y_axis_values = np.zeros(len(bins)-1)
        y_axis_error = np.zeros(len(bins)-1)

        for i2, bin in enumerate(list(bins[1:])):
            not_lower_than = local_x_values > bins[i2]
            not_higher_than = local_x_values < bins[i2+1]
            mask = not_higher_than & not_lower_than
            bin_ion_count = np.count_nonzero(mask)

            # bin_ion_std = np.std(local_comp[mask])
            # try:
                # bin_error = bin_ion_std / (len(local_comp[mask])**0.5) # SEM
            # except ZeroDivisionError:
                # bin_error = 0
            y_axis_values[i2] = bin_ion_count / df_comparison["Volume [nm3]"][i]
            # y_axis_error[i2] = bin_error

        fig.add_trace(go.Bar(
            x=bins,
            y=y_axis_values, 
            name=sample_titles[i], 
            marker_color=bar_colors[i],
            # mode="markers",
            # opacity=1,
            # marker_line_width=20,
            # error_y=dict(
            #     type='data', # value of error bar given in data coordinates
            #     # array=Cu_comp_std/(Cu_count**0.5), # error = SEM based on sigma / sqrt(N)
            #     # array=y_axis_error, # error = sigma
            #     visible=True,
            #     width=1
            #     )
        ))
        # break

    # Overlay both histograms
    fig.update_layout(
        yaxis_title=f"Number density [clusters / nm-3]",
        title_text='Size Distribution',
        barmode='group' # group together boxes of the different traces for each value of x
    )

    # set up the x axis
    ticktext = []
    tickvals = []
    for i in range(len(bins[1:])):
        text = "{:.1f}-{:.1f}".format(bins[i], bins[i+1])
        ticktext.append(text)
        tickvals.append(bins[i]) 

    fig.update_xaxes(
        ticktext=ticktext,
        tickvals=tickvals,
        title_text=x_axis
    )

    fig.show()
    return None


def makeAxis(title, tickangle):
    return {
    'title': title,
    'titlefont': { 'size': 15 },
    'tickangle': tickangle,
    'tickfont': { 'size': 15 },
    'tickcolor': 'rgba(0,0,0,0)',
    'ticklen': 5,
    'showline': True,
    'showgrid': True
    }


def plot_comparison_binned_ternary_binned(
    cluster_dataframes:dict_with_pandas_dfs,
    ternary_elements:string_vector,
    bins:float_vector,
    sample_titles:Optional[string_vector]=None
) -> None:

    # plot mean compositions on ternary diagrams
    fig = go.Figure()
    # set up colors for each sample
    bar_colors = chart_colors(len(cluster_dataframes))

    # ion = "Cu"
    # bins = np.arange(0, 2, 0.3) # ONE STEP LONGER THAN YOU WANT IT TO BE
    # ternary_elements = ["Mn+Si", "Ni", "Cu"] # ONLY THREE
    # check if there are only 3 elements
    if len(ternary_elements) != 3:
        print("You have to specify exactly three elements / columns!")

    # all_ternary_dfs = []

    for i in cluster_dataframes:
        
        # find ternary elements, and create sums if required
        for ternary_element in ternary_elements:
            if "+" in ternary_element:
                # create a list of elements and sum them up
                sub_ternary_elements = ternary_element.split(sep="+")
                cluster_dataframes[i][ternary_element] = cluster_dataframes[i][sub_ternary_elements].sum(axis=1)


        # find bin stats to smooth the graphs and generalise the data
        
        # local_cs_normalised = local_cs / etas[i] # cluster size / LEAP efficiency
        local_rgyration = cluster_dataframes[i]["r_gyration"]

        # create new columns if needed
        # dfs[i]["Mn+Si"] = dfs[i]["Si"] + dfs[i]["Mn"]

        # iterate through every bin and find mean and std of ion
        y_axis_values = [0, 0, 0]
        y_axis_values[0] = np.zeros(len(bins)-1)
        y_axis_values[1] = np.zeros(len(bins)-1)
        y_axis_values[2] = np.zeros(len(bins)-1)

        y_axis_error = [0, 0, 0]
        y_axis_error[0] = np.zeros(len(bins)-1)
        y_axis_error[1] = np.zeros(len(bins)-1)
        y_axis_error[2] = np.zeros(len(bins)-1)

        bin_cluster_count = np.zeros(len(bins)-1)

        for i2, bin in enumerate(list(bins[1:])):
            not_lower_than = local_rgyration > bins[i2]
            not_higher_than = local_rgyration < bins[i2+1]
            mask = not_higher_than & not_lower_than
        
            # iterate through three ions to plot
            for i3, ion in enumerate(ternary_elements):
                
                bin_ion_mean = np.mean(cluster_dataframes[i][ion][mask])
                bin_ion_std = np.std(cluster_dataframes[i][ion][mask])
                try:
                    bin_error = bin_ion_std / (len(cluster_dataframes[i][ion][mask])**0.5) # SEM
                except ZeroDivisionError:
                    bin_error = 0
                y_axis_values[i3][i2] = bin_ion_mean
                y_axis_error[i3][i2] = bin_error

            # find cluster count in each bin
            bin_cluster_count[i2] = len(cluster_dataframes[i][ion][mask])


        # move all the data to new pandas df, it's easier to plot this way
        this_df = pandas.DataFrame()
        
        this_df["Cluster Radius"] = bins[:-1]
        
        this_df[ternary_elements[0]] = y_axis_values[0]
        this_df[ternary_elements[1]] = y_axis_values[1]
        this_df[ternary_elements[2]] = y_axis_values[2]

        this_df["Bin Cluster Count"] = bin_cluster_count
        column_name = "{} Error"
        this_df[column_name.format(ternary_elements[0])] = y_axis_error[0]
        this_df[column_name.format(ternary_elements[1])] = y_axis_error[1]
        this_df[column_name.format(ternary_elements[2])] = y_axis_error[2]

        # add sample column
        this_df["Sample"] = np.array(([sample_titles[i]] * len(bin_cluster_count))) 

        # add informative cluster size ranges for hover id
        ticktext = []
        tickvals = []
        for i2 in range(len(bins[1:])):
            text = "{:.1f}-{:.1f}".format(bins[i2], bins[i2+1])
            ticktext.append(text)
            tickvals.append(bins[i2])
        this_df["Radius Range [nm]"] = ticktext

        # add to the list of all dataframes and concatenate later
        # all_ternary_dfs.append(this_df)

        # prepare marker sizes for these bins
        how_many_markers = len(this_df[ternary_elements[0]])
        marker_sizes = np.linspace(8, 14, num=how_many_markers)
        marker_sizes = np.round(marker_sizes)

        # plot average compositions for this sample
        fig.add_trace(go.Scatterternary(
            a=this_df[ternary_elements[0]], 
            b=this_df[ternary_elements[1]], 
            c=this_df[ternary_elements[2]], 
            opacity=0.7,
            mode="lines+markers",
            name=sample_titles[i],
            marker={
                # 'symbol': 100,
                'color': bar_colors[i],
                'size': marker_sizes # TODO improve! make it work with different lengths of bins variable
                }
        ))

    fig.update_layout({
        'ternary': {
            # 'sum': 100,
            'aaxis': makeAxis(ternary_elements[0], 0),
            'baxis': makeAxis(ternary_elements[1], 45),
            'caxis': makeAxis(ternary_elements[2], -45)
        },
        'annotations': [{
        'showarrow': False,
        'text': 'How concentration changes with size',
            'x': 0.5,
            'y': 1.3,
            'font': { 'size': 15 }
        }]
    })

    fig.show()

    return None


