## R ggplot theme

* Black font.
* No grid line, keep x and y axis.
* Adjust line width.

``` R
# element_blank() if not show
theme_classic() +
theme(
  legend.position = "right",
  plot.title = element_text(hjust = 0.5, color = "black"),
  axis.title = element_text(color = "black", size = 12),
  axis.text.y = element_text(color = "black", size = 12),
  axis.text.x = element_text(color = "black", size = 12),
  legend.text = element_text(color = "black", size = 12),
  legend.title = element_text(color = "black", size = 12),
  panel.grid = element_blank(),
  axis.line.x = element_line(color = "black", linewidth = 1, lineend = "square"),
  axis.line.y = element_line(color = "black", linewidth = 1, lineend = "square"),
  axis.ticks.x = element_line(color = "black", linewidth = 1),
  axis.ticks.y = element_line(color = "black", linewidth = 1),
  axis.ticks.length = unit(12, "pt")
)
```


## Python matplotlib theme
### plot header
Add this header, so that the font could be edited in Ai.
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
