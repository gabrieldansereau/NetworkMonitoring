"""
    heatmapcb(mat; label="", kw...)

Just a heatmap with colorbar with a colorbar for quick visualization. Because
doesn't it make sense?
"""
function heatmapcb(mat; label="", kw...)
    f, ax, p = heatmap(mat; kw...)
    Colorbar(f[1, end + 1], p; label=label)
    return f
end