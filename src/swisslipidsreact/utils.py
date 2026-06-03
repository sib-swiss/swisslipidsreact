import os

# Use environment variable SLR_DEBUG to set debug levels for logging and other things.
DEBUG = int( os.getenv( "SLR_DEBUG", 0 ) )

stats = {
    "rhea": [], # unique Rhea IDs
    "reactions": [], # enumerated reactions
    "participants": [],
    "lipids": []
}

def add_statistics(logger, name: str, message: str):

    stats[name].append( message )
    logger.info("Statistics: %s", message)

def save_statistics(output_dir):

    f_results_statistics = "statistics.txt"
    with open(os.path.join(output_dir, f_results_statistics), 'w') as f:
        
        f.write("\n\nRhea ID statistics:\n\n")
        for line in stats['rhea']:
            f.write(line)

        f.write("\n\nRhea participants statistics:\n\n")
        for line in stats['participants']:
            f.write(line)

        f.write("\n\nSLM lipids statistics:\n\n")
        for line in stats['lipids']:
            f.write(line)

        f.write("\n\nEnumerated reactions statistics:\n\n")
        for line in stats['reactions']:
            f.write(line)


def debug_df_first_row( df, separator='\n' ):
    """
    Shows the index, column name and value for all columns of a dataframe's first row.
    Use for debugging to understand what is stored in a dataframe.
    """
    columns = df.columns
    first_row = next( df.itertuples( index=False, name=None ) )
    row_string = separator.join( f'[{i}] {columns[i]}  ==>  {v}' for i, v in enumerate( first_row ) )
    return row_string
