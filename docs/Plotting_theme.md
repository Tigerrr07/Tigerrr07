## R ggplot theme

* Black font.
* No grid line, keep x and y axis.
* Adjust line width.
* Adjust ticks width.

``` R
# element_blank() if not show
theme_classic() +
theme(
  panel.grid = element_blank(),
  legend.position = "right",
  plot.title = element_text(hjust = 0.5, color = "black"),
  legend.text = element_text(color = "black", size = 14),
  legend.title = element_text(color = "black", size = 14),
  axis.title = element_text(color = "black", size = 14),
  axis.text.x = element_text(color = "black", size = 14),
  axis.text.y = element_text(color = "black", size = 14),
  axis.line.x = element_line(color = "black", linewidth = 0.7, lineend = "square"),
  axis.line.y = element_line(color = "black", linewidth = 0.7, lineend = "square"),
  axis.ticks.x = element_line(color = "black", linewidth = 0.7),
  axis.ticks.y = element_line(color = "black", linewidth = 0.7),
  axis.ticks.length = unit(4, "pt")
)
```


## Python matplotlib theme
### plot header
Add this header, so that the font could be edited in adobe illustrator.
``` python
from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['pdf.use14corefonts'] = True
rcParams['pdf.fonttype'] = 42 
rcParams['ps.fonttype'] = 42
```

### plot details
* No grid line, keep x and y axis.
* Adjust line width.

``` python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=figsize)

# Keep left and bottom line
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)

# Enlarge ticks width and length
ax.tick_params(axis='x', width=3, length=6)
ax.tick_params(axis='y', width=3, length=6)

# If not show ticks
ax.set_xticks([])
ax.set_yticks([])
fig.savefig(fig_name) # pdf or svg
```

### Plot scanpy spatial maps
* Save given size of spatial maps
``` python
# if ct is cell_type from deconvolution
with plt.rc_context({'font.size': 15, 'legend.markerscale': 2.0}):
    sc.pl.spatial(
        adata, 
        color=ct, 
        spot_size=spot_size, 
        frameon=False, 
        cmap='magma', 
        img_key=None,
        vmin=0, 
        vmax='p99', 
        show=False, 
        title='',
        colorbar_loc=None
    )
    
    # Get current graph
    fig = plt.gcf()    
    fig.set_size_inches(4, 4)
    fig.subplots_adjust(left=0.001, right=0.999, top=0.999, bottom=0.001) # adjust
    plt.savefig(f"{ct.replace(' ', '_')}.pdf")
    plt.close()
    
    # Create separate colobar
    mappable = fig.axes[0].collections[0]
    fig_cbar, ax_cbar = plt.subplots(figsize=(0.6, 2))
    fig_cbar.subplots_adjust(left=0.1, right=0.3, top=0.95, bottom=0.05)
    cbar = plt.colorbar(mappable, cax=ax_cbar)
    cbar.ax.tick_params(labelsize=15)
    
    # save colorbar
    plt.savefig(f"{ct.replace(' ', '_')}_colorbar.pdf")
    plt.close(fig_cbar)
```


## Change defult font in Python
Copy system Arial.tff to a remote path:

``` bash
~/.local/share/fonts/arial
```

Change the settings
``` python
import matplotlib.pyplot as plt
from matplotlib import font_manager

font_dir = "~/.local/share/fonts/arial"
for f in font_manager.findSystemFonts(fontpaths=[font_dir]):
    font_manager.fontManager.addfont(f)

plt.rcParams["font.family"] = "Arial"
plt.rcParams["axes.unicode_minus"] = False


# fc-cache -f -v ~/.local/share/fonts
# fc-match Arial
```