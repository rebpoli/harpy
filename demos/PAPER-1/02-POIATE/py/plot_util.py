
import os

#
# Save write to avoid NFS issues
#
def savefig( fig, fn ) :
    print(f"Saving {fn} ...")
    import shutil
    bn = os.path.basename(fn)
    temp_file = f'/tmp/{bn}'
    fig.savefig(temp_file, dpi=300)
    shutil.move(temp_file, fn)  # Atomic operation    
