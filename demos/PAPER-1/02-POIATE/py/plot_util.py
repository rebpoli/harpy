
#
# Save write to avoid NFS issues
#
def savefig( fig, fn ) :
    print(f"Saving {fn} ...")
    import shutil
    temp_file = f'/tmp/{fn}.png'
    fig.savefig(temp_file, dpi=300)
    shutil.move(temp_file, fn)  # Atomic operation    
