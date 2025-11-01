import subprocess
import os
from pathlib import Path
import shutil


# Path to your METISSE executable and working directory

METISSE_EXE = "./metisse"

def compile_metisse(METISSE_DIR,run_dir):
    """
    Compile METISSE by running the local 'mk' build script in METISSE_DIR,
    copy the produced 'metisse' executable into run_dir, then run './clean'.

    Parameters
    ----------
    METISSE_DIR : str or Path
        Directory containing the METISSE source and the 'mk' build script.
    run_dir : str or Path
        Directory where the compiled 'metisse' executable will be copied.

    """
    
    # Run the mk build script using the shell. cwd sets the working directory.
    result = subprocess.run("./mk",
                            cwd=METISSE_DIR,
                            stdout=subprocess.PIPE,
                            shell=True, text=True)
    
    # If the build failed (non-zero exit), raise with stderr text.
    if result.returncode != 0:
        raise RuntimeError(f"METISSE compilation failed:\n{result.stderr}")

    # Copy the built executable named "metisse" from METISSE_DIR into run_dir
    # shutil.copy2 preserves metadata (timestamps, permissions when possible).
    shutil.copy2(os.path.join(METISSE_DIR, "metisse"), run_dir)

    # Run a clean step in the source directory 
    subprocess.run("./clean",
                            cwd=METISSE_DIR,
                            stdout=subprocess.PIPE,
                            shell=True, text=True)
    
    return 


def inlist_content(section_name,params):
    """Construct a Fortran-style namelist block (as text) from a dict.

    Parameters
    ----------
    section_name : str
        The namelist header, e.g. '&METISSE_input_controls'
    params : dict
        Mapping of parameter name -> Python value. Supported value types:
        - bool -> converted to .true. / .false.
        - str  -> wrapped in single quotes
        - numeric -> converted to string as-is

    Returns
    -------
    str
        Text of the namelist block, terminated with '/' and a trailing blank line.
    """

    lines = [section_name]
    for key, val in params.items():
        if isinstance(val, bool):
            # Fortran expects .true./.false. for booleans
            val_str = ".true." if val else ".false."
        elif isinstance(val, str):
            # Quote strings for Fortran
            val_str = f"'{val}'"
        else:
            # Numbers (ints/floats) are converted to their Python string repr
            val_str = str(val)
        lines.append(f"    {key} = {val_str}")
    # End the namelist block
    lines.append("/")
    # extra blank line, without it inlists give error
    lines.append(" ")
    text = "\n".join(lines)

    return text

def run_metisse(run_dir,main_params,metisse_params):
    """
    Launch the METISSE executable in run_dir and stream its stdout to the notebook.

    Parameters
    ----------
    run_dir : str or Path
        Directory where the metisse executable and inlists live; this becomes the cwd.
    main_params : dict
            Parameters for the '&SSE_input_controls' block.
    metisse_params : dict
            Parameters for the '&METISSE_input_controls' block.

    Raises
    ------
    RuntimeError
        If the METISSE process exits with a non-zero return code.
    """
    
    # Pass the params and write inlists
    # Compose and write the METISSE.input namelist
    path = Path(os.path.join(run_dir, "metisse.input"))
    text = inlist_content('&METISSE_input_controls',metisse_params)
    path.write_text(text)

    # Compose and write the main.input namelist
    path = Path(os.path.join(run_dir, "main.input"))
    text = inlist_content('&SSE_input_controls',main_params)
    path.write_text(text)
                 
    proc = subprocess.Popen([(METISSE_EXE)],
                            cwd=str(run_dir),
                            shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True)
    # Stream stdout line-by-line.
    for line in iter(proc.stdout.readline, ''):
        print(line, end='')

    # Wait for the process to finish
    proc.wait()

    # On non-zero exit, read stderr and raise
    if proc.returncode != 0:
        err = proc.stderr.read()
        raise RuntimeError(f"METISSE run failed:\n{err}")
    
    return 
