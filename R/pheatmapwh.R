## The code in this file is based on that in the file pheatmap.r in the CRAN package pheatmap by  Raivo Kolde <rkolde@gmail.com>,
##    with minor (but important to us) formatting modifications.

## These global parameters were introduced by Wolfgang Huber. In Raivo's code, they were hardcoded numbers
gaps_width <- c( col = 26, row = 10)
horizonal_annot_mat <- 0 # 8

lo = function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, treeheight_col, treeheight_row, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, gaps_row, gaps_col, ...){
    # Get height of colnames and length of rownames
    if(!is.null(coln[1])){
        t = c(coln, colnames(annotation_row))
        longest_coln = which.max(strwidth(t, units = 'in'))
        gp = list(fontsize = fontsize_col, ...)
        coln_height = unit(1, "grobheight", textGrob(t[longest_coln], rot = 90, gp = do.call(gpar, gp))) + unit(10, "bigpts")
    }
    else{
        coln_height = unit(5, "bigpts")
    }
    
    if(!is.null(rown[1])){
        t = c(rown, colnames(annotation_col))
        longest_rown = which.max(strwidth(t, units = 'in'))
        gp = list(fontsize = fontsize_row, ...)
        rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], gp = do.call(gpar, gp))) ## + unit(10, "bigpts")
    }
    else{
        rown_width = unit(5, "bigpts")
    }
    
    gp = list(fontsize = fontsize, ...)
    # Legend position
    if(!is.na2(legend)){
        longest_break = which.max(nchar(names(legend)))
        longest_break = unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
        title_length = unit(1.1, "grobwidth", textGrob("Scale", gp = gpar(fontface = "bold", ...)))
        legend_width = unit(12, "bigpts") + longest_break * 1.2
        legend_width = max(title_length, legend_width)
    }
    else{
        legend_width = unit(0, "bigpts")
    }
    
    # Set main title height
    if(is.na(main)){
        main_height = unit(0, "npc")
    }
    else{
        main_height = unit(1.5, "grobheight", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }
    
    # Column annotations
    textheight = unit(fontsize, "bigpts")
    
    if(!is.na2(annotation_col)){
        # Column annotation height 
        annot_col_height = ncol(annotation_col) * (textheight + unit(2, "bigpts")) + unit(horizonal_annot_mat, "bigpts")
        
        # Width of the correponding legend
        t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col)) 
        annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
        if(!annotation_legend){
            annot_col_legend_width = unit(0, "npc")
        }
    }
    else{
        annot_col_height = unit(0, "bigpts")
        annot_col_legend_width = unit(0, "bigpts")
    }
    
    # Row annotations
    if(!is.na2(annotation_row)){
        # Row annotation width 
        annot_row_width = ncol(annotation_row) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")
        
        # Width of the correponding legend
        t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row)) 
        annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
        if(!annotation_legend){
            annot_row_legend_width = unit(0, "npc")
        }
    }
    else{
        annot_row_width = unit(0, "bigpts")
        annot_row_legend_width = unit(0, "bigpts")
    }
    
    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
    
    # Tree height
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts") 
    
    # Set cell sizes
    mat_width = if(is.na(cellwidth)){
        unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_row_width - annot_legend_width 
    } else{
        unit(cellwidth * ncol, "bigpts") + length(gaps_col) * unit(gaps_width["col"], "bigpts")
    }
    
    mat_height = if(is.na(cellheight)){
        unit(1, "npc") - main_height - coln_height - treeheight_col - annot_col_height
    } else{
        unit(cellheight * nrow, "bigpts") + length(gaps_row) * unit(gaps_width["row"], "bigpts")
    }    
    
    # Produce gtable
    gt = gtable(widths = unit.c(             treeheight_row, annot_row_width,   mat_width, rown_width, legend_width, annot_legend_width),
               heights = unit.c(main_height, treeheight_col, annot_col_height, mat_height, coln_height), vp = viewport(gp = do.call(gpar, gp)))
    
    cw = convertWidth( mat_width  - (length(gaps_col) * unit(gaps_width["col"], "bigpts")), "bigpts", valueOnly = TRUE) / ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(gaps_width["row"], "bigpts")), "bigpts", valueOnly = TRUE) / nrow
    
    # Return minimal cell dimension in bigpts to decide if borders are drawn
    mindim = min(cw, ch) 
    
    res = list(gt = gt, mindim = mindim)
    
    return(res)
}

find_coordinates = function(n, gaps, m = 1:n, what) {
    if (length(gaps) == 0)
        return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc") ))

    if(max(gaps) > n)
        stop("Gaps do not match with matrix size")

    size = (1 / n) * (unit(1, "npc") - length(gaps) * unit(gaps_width[what], "bigpts"))
    
    gaps2 = apply(sapply(gaps, function(gap, x){x > gap}, m), 1, sum) 
    coord = m * size + (gaps2 * unit(gaps_width[what], "bigpts"))
    
    return(list(coord = coord, size = size))
}

draw_dendrogram = function(hc, gaps, horizontal = TRUE){
    h = hc$height / max(hc$height) / 1.05
    m = hc$merge
    o = hc$order
    n = length(o)

    m[m > 0] = n + m[m > 0] 
    m[m < 0] = abs(m[m < 0])

    dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
    dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)

    for(i in 1:nrow(m)){
        dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
        dist[n + i, 2] = h[i]
    }
    
    draw_connection = function(x1, x2, y1, y2, y){
        res = list(
            x = c(x1, x1, x2, x2),
            y = c(y1, y, y, y2)
        )
        
        return(res)
    }
    
    x = rep(NA, nrow(m) * 4)
    y = rep(NA, nrow(m) * 4)
    id = rep(1:nrow(m), rep(4, nrow(m)))
    
    for(i in 1:nrow(m)){
        c = draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
        k = (i - 1) * 4 + 1
        x[k : (k + 3)] = c$x
        y[k : (k + 3)] = c$y
    }
    
    x = find_coordinates(n, gaps, x * n, what = ifelse(horizontal, "col", "row"))$coord
    y = unit(y, "npc")
    
    if(!horizontal){
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a
    }
    res = polylineGrob(x = x, y = y, id = id)
    
    return(res)
}

draw_matrix = function(matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color){
    n = nrow(matrix)
    m = ncol(matrix)
    
    coord_x = find_coordinates(m, gaps_cols, what = "col")
    coord_y = find_coordinates(n, gaps_rows, what = "row")
    
    x = coord_x$coord - 0.5 * coord_x$size
    y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)
    
    coord = expand.grid(y = y, x = x)
    
    res = gList()
    
    res[["rect"]] = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = coord_y$size, gp = gpar(fill = matrix, col = border_color))
    
    if(attr(fmat, "draw")){
        res[["text"]] = textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
    }
    
    res = gTree(children = res)
    
    return(res)
}

draw_colnames = function(coln, gaps, ...){
    coord = find_coordinates(length(coln), gaps, what = "col")
    x = coord$coord - 0.5 * coord$size
    
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))
    
    return(res)
}

draw_rownames = function(rown, gaps, ...){
    coord = find_coordinates(length(rown), gaps, what = "row")
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
    
    res = textGrob(rown, x = unit(3, "bigpts"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))
    
    return(res)
}

draw_legend = function(color, breaks, legend, ...){
    height = min(unit(1, "npc"), unit(150, "bigpts"))
    
    legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
    legend_pos = height * legend_pos + (unit(1, "npc") - height)
    
    breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
    breaks = height * breaks + (unit(1, "npc") - height)
    
    h = breaks[-1] - breaks[-length(breaks)]
    
    rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
    text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
    
    res = grobTree(rect, text)
    
    return(res)
}

convert_annotations = function(annotation, annotation_colors){
    new = annotation
    for(i in seq_len(ncol(annotation))) {
        a = annotation[, i]
        b = annotation_colors[[colnames(annotation)[i]]]
        if (is.character(a) || is.factor(a)){
            a = as.character(a)
            
            if(length(setdiff(setdiff(a, NA), names(b))) > 0){
                stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
            }
            new[, i] = b[a]
        }
        else{
            a = cut(a, breaks = 100)
            new[, i] = colorRampPalette(b)(100)[a]
        }
    }
    return(as.matrix(new))
}

draw_annotations = function(converted_annotations, border_color, gaps, fontsize, horizontal){
    n = ncol(converted_annotations)
    m = nrow(converted_annotations)
    
    coord_x = find_coordinates(m, gaps, what = "col")
    
    x = coord_x$coord - 0.5 * coord_x$size
    
    # y = cumsum(rep(fontsize, n)) - 4 + cumsum(rep(2, n))
    y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1 
    y = unit(y, "bigpts")
    
    if (horizontal) {
        coord = expand.grid(x = x, y = y)
        res = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = unit(fontsize, "bigpts"), gp = gpar(fill = converted_annotations, col = border_color))
    } else{
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a
        
        coord = expand.grid(y = y, x = x)
        res = rectGrob(x = coord$x, y = coord$y, width = unit(fontsize, "bigpts"), height = coord_x$size, gp = gpar(fill = converted_annotations, col = border_color))
    }
    
    return(res)
}

draw_annotation_names = function(annotations, fontsize, horizontal){
    n = ncol(annotations)
    
    x = unit(3, "bigpts")
    
    y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1 
    y = unit(y, "bigpts")
    
    if(horizontal){
        res = textGrob(colnames(annotations), x = x, y = y, hjust = 0, gp = gpar(fontsize = fontsize, fontface = 2))
    }
    else{
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a
        
        res = textGrob(colnames(annotations), x = x, y = y, vjust = 0.5, hjust = 0, rot = 270, gp = gpar(fontsize = fontsize, fontface = 2))
    }
    
    return(res)
}

draw_annotation_legend = function(annotation, annotation_colors, border_color, ...){
    y = unit(1, "npc")
    text_height = unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))
    
    res = gList()

    ## by WH: remove those with leading blank, or "n.d."
    columnsToDraw <- names(annotation)
    columnsToDraw <- columnsToDraw[ !grepl("^ ", columnsToDraw) ]
    
    for(i in columnsToDraw) {
        res[[i]] = textGrob(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = "bold", ...))
        
        y = y - 1.5 * text_height
        if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
            aci <- annotation_colors[[i]]
            aci <- aci[ names(aci) != "n.d." ]
            n <- length(aci)
            yy <- y - (seq_len(n) - 1) * 2 * text_height
            
            res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = 2 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = aci))
            res[[paste(i, "t")]] = textGrob(names(aci), x = text_height * 2.4, y = yy - text_height, hjust = 0, vjust = 0.5, gp = gpar(...))
            
            y = y - n * 2 * text_height
            
        }
        else{
            yy = y - 8 * text_height + seq(0, 1, 0.25)[-1] * 8 * text_height
            h = 8 * text_height * 0.25
            
            res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = 2 * text_height, gp = gpar(col = NA, fill = colorRampPalette(annotation_colors[[i]])(4)))
            res[[paste(i, "r2")]] = rectGrob(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = 8 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = NA))
            
            txt = rev(range(grid.pretty(range(annotation[[i]], na.rm = TRUE))))
            yy = y - c(1, 7) * text_height
            res[[paste(i, "t")]]  = textGrob(txt, x = text_height * 2.4, y = yy, hjust = 0, vjust = 0.5, gp = gpar(...))
            y = y - 8 * text_height
        }
        y = y - 1.5 * text_height
    }
    
    res = gTree(children = res)
    
    return(res)
}

draw_main = function(text, ...){
    res = textGrob(text, gp = gpar(fontface = "bold", ...))
    
    return(res)
}

vplayout = function(x, y){
    return(viewport(layout.pos.row = x, layout.pos.col = y))
}

heatmap_motor = function(matrix, border_color, cellwidth, cellheight, tree_col, tree_row, treeheight_col, treeheight_row, filename, width, height, breaks, color, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, fmat, fontsize_number, number_color, gaps_col, gaps_row, labels_row, labels_col, ...){
    # Set layout
    lo = lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, gaps_row = gaps_row, gaps_col = gaps_col,  ...)
    
    res = lo$gt
    mindim = lo$mindim
    
    if(!is.na(filename)){
        if(is.na(height)){
            height = convertHeight(gtable_height(res), "inches", valueOnly = T)
        }
        if(is.na(width)){
            width = convertWidth(gtable_width(res), "inches", valueOnly = T)
        }
        
        # Get file type
        r = regexpr("\\.[a-zA-Z]*$", filename)
        if(r == -1) stop("Improper filename")
        ending = substr(filename, r + 1, r + attr(r, "match.length"))

        f = switch(ending,
            pdf = function(x, ...) pdf(x, ...),
            png = function(x, ...) png(x, units = "in", res = 300, ...),
            jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
            jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
            tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
            bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
            stop("File type should be: pdf, png, bmp, jpg, tiff")
        )
        
        # print(sprintf("height:%f width:%f", height, width))
        
        # gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)

        f(filename, height = height, width = width)
        gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)
        grid.draw(gt)
        dev.off()
        
        return(gt)
    }
    
    # Omit border color if cell size is too small 
    if(mindim < 3) border_color = NA
    
    # Draw title
    if(!is.na(main)){
        elem = draw_main(main, fontsize = 1.3 * fontsize, ...)
        res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main")
    }
    
    # Draw tree for the columns
    if(!is.na2(tree_col) & treeheight_col != 0){
        elem = draw_dendrogram(tree_col, gaps_col, horizontal = T)
        res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
    }
    
    # Draw tree for the rows
    if(!is.na2(tree_row) & treeheight_row != 0){
        elem = draw_dendrogram(tree_row, gaps_row, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }
    
    # Draw matrix
    elem = draw_matrix(matrix, border_color, gaps_row, gaps_col, fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", name = "matrix")
    
    # Draw colnames
    if(length(labels_col) != 0){
        pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, ...)
        elem = do.call(draw_colnames, pars)
        res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", name = "col_names")
    }
    
    # Draw rownames
    if(length(labels_row) != 0){
        pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, ...)
        elem = do.call(draw_rownames, pars)
        res = gtable_add_grob(res, elem, t = 4, l = 4, clip = "off", name = "row_names")
    }
    
    # Draw annotation tracks on cols
    if(!is.na2(annotation_col)){
        # Draw tracks
        converted_annotation = convert_annotations(annotation_col, annotation_colors)
        elem = draw_annotations(converted_annotation, border_color, gaps_col, fontsize, horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", name = "col_annotation")
        
        # Draw names
        elem = draw_annotation_names(annotation_col, fontsize, horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", name = "row_annotation_names")
        
    }
    
    # Draw annotation tracks on rows
    if(!is.na2(annotation_row)){
        # Draw tracks
        converted_annotation = convert_annotations(annotation_row, annotation_colors)
        elem = draw_annotations(converted_annotation, border_color, gaps_row, fontsize, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", name = "row_annotation")
        
        # Draw names
        if(length(labels_col) != 0){
            elem = draw_annotation_names(annotation_row, fontsize, horizontal = F)
            res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", name = "row_annotation_names")
        }
    }
    
    # Draw annotation legend
    annotation = c(annotation_col[length(annotation_col):1], annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]
    
    if(length(annotation) > 0 & annotation_legend){
        elem = draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)
        
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, clip = "off", name = "annotation_legend")
    }
    
    # Draw legend
    if(!is.na2(legend)){
        elem = draw_legend(color, breaks, legend, fontsize = fontsize, ...)
        
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, clip = "off", name = "legend")
    }
    
    return(res)
}

generate_breaks = function(x, n, center = F){
    if(center){
        m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
        res = seq(-m, m, length.out = n + 1)
    }
    else{
        res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
    }
    
    return(res)
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
    return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

scale_colours = function(mat, col = rainbow(10), breaks = NA){
    mat = as.matrix(mat)
    return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

cluster_mat = function(mat, distance, method){
    if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
        stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
    }
    if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
        stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
    }
    if(distance[1] == "correlation"){
        d = as.dist(1 - cor(t(mat)))
    }
    else{
        if(class(distance) == "dist"){
            d = distance
        }
        else{
            d = dist(mat, method = distance)
        }
    }
    
    return(hclust(d, method = method))
}

scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

scale_mat = function(mat, scale){
    if(!(scale %in% c("none", "row", "column"))){
        stop("scale argument shoud take values: 'none', 'row' or 'column'")
    }
    mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
    return(mat)
}

generate_annotation_colours = function(annotation, annotation_colors, drop){
    if(is.na2(annotation_colors)){
        annotation_colors = list()
    }
    count = 0
    for(i in 1:length(annotation)){
        annotation[[i]] = annotation[[i]][!is.na(annotation[[i]])]
        if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
            if (is.factor(annotation[[i]]) & !drop){
                count = count + length(levels(annotation[[i]]))
            }
            else{
                count = count + length(unique(annotation[[i]]))
            }
        }
    }
    
    factor_colors = dscale(factor(1:count), hue_pal(l = 75))
    
    set.seed(3453)
    
    cont_counter = 2
    for(i in 1:length(annotation)){
        if(!(names(annotation)[i] %in% names(annotation_colors))){
            if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
                n = length(unique(annotation[[i]]))
                if (is.factor(annotation[[i]]) & !drop){
                    n = length(levels(annotation[[i]]))
                }
                ind = sample(1:length(factor_colors), n)
                annotation_colors[[names(annotation)[i]]] = factor_colors[ind]
                l = levels(as.factor(annotation[[i]]))
                l = l[l %in% unique(annotation[[i]])]
                if (is.factor(annotation[[i]]) & !drop){
                    l = levels(annotation[[i]])
                }
                
                names(annotation_colors[[names(annotation)[i]]]) = l
                factor_colors = factor_colors[-ind]
            }
            else{
                annotation_colors[[names(annotation)[i]]] = brewer_pal("seq", cont_counter)(5)[1:4]
                cont_counter = cont_counter + 1
            }
        }
    }
    return(annotation_colors)
}

kmeans_pheatmap = function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...){
    # Filter data
    if(!is.na(sd_limit)){
        s = apply(mat, 1, sd)
        mat = mat[s > sd_limit, ]    
    }
    
    # Cluster data
    set.seed(1245678)
    km = kmeans(mat, k, iter.max = 100)
    mat2 = km$centers
    
    # Compose rownames
    t = table(km$cluster)
    rownames(mat2) = sprintf("cl%s_size_%d", names(t), t)
    
    # Draw heatmap
    pheatmapwh(mat2, ...)
}

find_gaps = function(tree, cutree_n){
    v = cutree(tree, cutree_n)[tree$order]
    gaps = which((v[-1] - v[-length(v)]) != 0)
    
}

is.na2 = function(x){
    if(is.list(x) | length(x) > 1){
        return(FALSE)
    }
    if(length(x) == 0){
        return(TRUE)
    }
    
    return(is.na(x))
}

identity2 = function(x, ...){
    return(x)
}

pheatmapwh = function(mat, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, border_color = "grey60", cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", clustering_callback = identity2, cutree_rows = NA, cutree_cols = NA,  treeheight_row = ifelse(cluster_rows, 50, 0), treeheight_col = ifelse(cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA, legend_labels = NA, annotation_row = NA, annotation_col = NA, annotation = NA, annotation_colors = NA, annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8 * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL, labels_col = NULL, filename = NA, width = NA, height = NA, silent = FALSE, ...){
    
    # Set labels
    if(is.null(labels_row)){
        labels_row = rownames(mat)
    }
    if(is.null(labels_col)){
        labels_col = colnames(mat)
    }
    
    # Preprocess matrix
    mat = as.matrix(mat)
    if(scale != "none"){
        mat = scale_mat(mat, scale)
        if(is.na2(breaks)){
            breaks = generate_breaks(mat, length(color), center = T)
        }
    }
    
    
    # Kmeans
    if(!is.na(kmeans_k)){
        # Cluster data
        km = kmeans(mat, kmeans_k, iter.max = 100)
        mat = km$centers
        
        # Compose rownames
        t = table(km$cluster)
        labels_row = sprintf("Cluster: %s Size: %d", names(t), t)
    }
    else{
        km = NA
    }
    
    # Format numbers to be displayed in cells
    if(is.matrix(display_numbers) | is.data.frame(display_numbers)){
        if(nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)){
            stop("If display_numbers provided as matrix, its dimensions have to match with mat")
        }
        
        display_numbers = as.matrix(display_numbers)
        fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
        fmat_draw = TRUE
        
    }
    else{
        if(display_numbers){
            fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
            fmat_draw = TRUE
        }
        else{
            fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
            fmat_draw = FALSE
        }
    }
    
    # Do clustering
    if(cluster_rows){
        tree_row = cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)
        tree_row = clustering_callback(tree_row, mat)
        mat = mat[tree_row$order, , drop = FALSE]
        fmat = fmat[tree_row$order, , drop = FALSE]
        labels_row = labels_row[tree_row$order]
        if(!is.na(cutree_rows)){
            gaps_row = find_gaps(tree_row, cutree_rows)
        }
        else{
            gaps_row = NULL
        }
    }
    else{
        tree_row = NA
        treeheight_row = 0
    }
    
    if(cluster_cols){
        tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)
        tree_col = clustering_callback(tree_col, t(mat))
        mat = mat[, tree_col$order, drop = FALSE]
        fmat = fmat[, tree_col$order, drop = FALSE]
        labels_col = labels_col[tree_col$order]
        if(!is.na(cutree_cols)){
            gaps_col = find_gaps(tree_col, cutree_cols)
        }
        else{
            gaps_col = NULL
        }
    }
    else{
        tree_col = NA
        treeheight_col = 0
    }
    
    attr(fmat, "draw") = fmat_draw
    
    # Colors and scales
    if(!is.na2(legend_breaks) & !is.na2(legend_labels)){
        if(length(legend_breaks) != length(legend_labels)){
            stop("Lengths of legend_breaks and legend_labels must be the same")
        }
    }
    
    
    if(is.na2(breaks)){
        breaks = generate_breaks(as.vector(mat), length(color))
    }
    if (legend & is.na2(legend_breaks)) {
        legend = grid.pretty(range(as.vector(breaks)))
        names(legend) = legend
    }
    else if(legend & !is.na2(legend_breaks)){
        legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
        
        if(!is.na2(legend_labels)){
            legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
            names(legend) = legend_labels
        }
        else{
            names(legend) = legend
        }
    }
    else {
        legend = NA
    }
    mat = scale_colours(mat, col = color, breaks = breaks)
    
    # Preparing annotations
    if(is.na2(annotation_col) & !is.na2(annotation)){
        annotation_col = annotation
    }
    # Select only the ones present in the matrix
    if(!is.na2(annotation_col)){
        annotation_col = annotation_col[colnames(mat), , drop = F]
    }
    
    if(!is.na2(annotation_row)){
        annotation_row = annotation_row[rownames(mat), , drop = F]
    }
    
    annotation = c(annotation_row, annotation_col)
    annotation = annotation[unlist(lapply(annotation, function(x) !is.na2(x)))]
    if(length(annotation) != 0){
        annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
    }
    else{
        annotation_colors = NA
    }
    
    if(!show_rownames){
        labels_row = NULL
    }
    
    if(!show_colnames){
        labels_col = NULL
    }
    
    # Draw heatmap
    gt = heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend, annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, number_color = number_color, gaps_row = gaps_row, gaps_col = gaps_col, labels_row = labels_row, labels_col = labels_col, ...)
    
    if(is.na(filename) & !silent){
        grid.newpage()
        grid.draw(gt)
    }
    
    invisible(list(tree_row = tree_row, tree_col = tree_col, kmeans = km, gtable = gt))
}


