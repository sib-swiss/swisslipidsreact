import os

# Use environment variable SLR_DEBUG to set debug levels for logging and other things.
DEBUG = int( os.getenv( "SLR_DEBUG", 0 ) )


def debug_df_first_row( df, separator='\n' ):
    """
    Shows the index, column name and value for all columns of a dataframe's first row.
    Use for debugging to understand what is stored in a dataframe.
    """
    columns = df.columns
    first_row = next( df.itertuples( index=False, name=None ) )
    row_string = separator.join( f'[{i}] {columns[i]}  ==>  {v}' for i, v in enumerate( first_row ) )
    return row_string
