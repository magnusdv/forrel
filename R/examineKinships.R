#' Check pedigree for relationship errors
#'
#' This function provides a convenient way to check for pedigree errors in a
#' linkage project or other situations where marker data is available for several
#' members. The function calls [IBDestimate()] to estimate IBD coefficients
#' for all indicated pairs of pedigree members and produces a colour-coded plot where wrong
#' relationships are easy to spot.
#'
#' @param x A `ped` object, or a list of such.
#' @param who A character vector of one or more of the words 'parents',
#'  'siblings', 'hugs' (= halfsibs/uncles/grandparents), 'cousins' and 'unrelated'.
#'  Two additional single-word values are possible: 'all' (all of the above, plus 'other')
#'  and 'close' (= 'parents', 'siblings', 'hugs', 'cousins').
#' @param interfam A character; either 'founders', 'none' or 'all', indicating
#'  which interfamiliar pairs of individuals should be included. Only relevant
#'  if `x` is a list of several `ped` objects.
#' @param makeplot A logical.
#' @param pch Plotting symbol (default: cross).
#' @param \dots Further plot arguments passed on to [showInTriangle()].
#'
#' @return A list of `data.frame`s (one for each relation category) with IBD estimates.
#' @author Magnus Dehli Vigeland
#' @seealso [IBDestimate()], [IBDtriangle()], [showInTriangle()]
#'
#' @examples
#' # Simulate 500 SNPs in a cousin pedigree
#' x = cousinPed(1)
#' x = simpleSim(x, 500, alleles = 1:2)
#' examineKinships(x)
#'
#' # Pretend we didn't know the brothers (3 and 6) were related
#' x1 = branch(x, 3)
#' x2 = branch(x, 6)
#' famid(x2) = 2
#'
#' # Notice the error: An 'unrelated' dot close to the sibling point
#' examineKinships(list(x1, x2))
#'
#'
examineKinships = function(x, who = "all", interfam = c("founders", "none", "all"),
                           makeplot = TRUE, pch = 4, ...) {

    pedList = is.pedList(x)

    groups = c("parents", "siblings", "hugs", "cousins", "other", "unrelated")
    group_txt = c("Parent-offspring", "Siblings", "Halfsib/Uncle/Grand", "1st cousins", "Other", "Unrelated")
    group_col = c("red", "blue", "green3", "magenta", "orange", "brown")
    names(group_txt) = names(group_col) = groups

    if (identical(who, "all")) {
        use_groups = groups
        if (pedList)
            use_groups = groups[groups != "other"]
    }
    else if (identical(who, "close"))
        use_groups = groups[1:4]
    else
        use_groups = sapply(who, match.arg, groups[-6])

    if (makeplot) {
        if (is.ped(x)) {
            g = .generations(x)
            marg = max(1, 12 - 2 * g)
            layout(rbind(1:2), widths = c(0.4, 0.6))
            plot(x, av = TRUE, cex = 0.8, title = "", margins = c(marg, 1, marg, 1))
        }
        else if (is.pedList(x))
            warnings("Automatic plotting is not implemented when 'x' is a list of ped objects. Use plotPedList().")

        IBDtriangle(relationships = c("UN", "PO", "MZ", "S", "H,U,G", "FC", "SC"))
        legend("topright", title = " According to pedigree:", title.adj = 0, pch = pch,
               lwd = 2, lty = NA, legend = group_txt[use_groups], col = group_col[use_groups])
    }

    P = CO = S = G = U = D = NULL
    p = s = g = co = d = u = NULL

    if ("parents" %in% use_groups) {
        P = related.pairs(x, "parents", available = TRUE)
        p = IBDestimate(x, P)
        showInTriangle(p$k0, p$k2, new = FALSE, col = group_col["parents"], pch = pch, ...)
    }
    if ("siblings" %in% use_groups) {
        S = rbind(related.pairs(x, "siblings", half = FALSE, available = TRUE))
        s = IBDestimate(x, S)
        showInTriangle(s$k0, s$k2, new = FALSE, col = group_col["siblings"], pch = pch, ...)
    }
    if ("cousins" %in% use_groups) {
        CO = rbind(related.pairs(x, "cousins", half = FALSE, available = TRUE),
                   related.pairs(x, "grandparents", degree = 3, available = TRUE),
                   related.pairs(x, "nephews_nieces", half = TRUE, available = TRUE),
                   related.pairs(x, "nephews_nieces", removal = 2, half = FALSE, available = TRUE))
        co = IBDestimate(x, CO)
        showInTriangle(co$k0, co$k2, new = FALSE, col = group_col["cousins"], pch = pch, ...)
    }
    if ("hugs" %in% use_groups) {
        G = rbind(related.pairs(x, "siblings", half = TRUE, available = TRUE),
                  related.pairs(x, "grandparents", available = TRUE),
                  related.pairs(x, "nephews_nieces", half = FALSE, available = TRUE))
        g = IBDestimate(x, G)
        showInTriangle(g$k0, g$k2, new = FALSE, col = group_col["hugs"], pch = pch, ...)
    }
    if ("unrelated" %in% use_groups) {
        interfam = match.arg(interfam)
        U = rbind(related.pairs(x, "unrelated", available = TRUE, interfam = interfam))
        u = IBDestimate(x, U)
        showInTriangle(u$k0, u$k2, new = FALSE, col = group_col["unrelated"], pch = pch, ...)
    }
    if ("other" %in% use_groups) {
        rest = t(combn(x$available, 2))
        taken = rbind(P, U, S, G, CO)
        D = rest[!(rest[, 1] * 1000 + rest[, 2]) %in% (taken[, 1] * 1000 + taken[, 2]), , drop = FALSE]
        d = IBDestimate(x, D)
        showInTriangle(d$k0, d$k2, new = FALSE, col = group_col["other"], pch = pch, ...)
    }
    res = list(parents = p, siblings = s, halfsibs_uncles_grandparents = g,
        cousins = co, other = d, unrelated = u)
    res = res[!sapply(res, is.null)]
    res
}

.IBDsuspects = function(x, radius) {
    # x output from examineKinships, radius = tolerated distance

    extract_distant = function(df, x0, y0, rad) {
        if (is.null(df))
            return(NULL)
        subs = subset(df, sqrt((df$k0 - x0)^2 + (df$k2 - y0)^2) > rad)
        if (nrow(subs) > 0)
            subs else NULL
    }

    susp = list(`suspects-parents` = extract_distant(x$parents, 0, 0, radius),
        `suspects-siblings` = extract_distant(x$siblings, 0.25, 0.25, radius),
        `suspects-halfsibs/grandparents/uncles` = extract_distant(x$halfsibs, 0.5, 0, radius),
        `suspects-cousins/halfuncles/granduncles` = extract_distant(x$cousins, 0.75, 0, radius),
        `suspects-other` = extract_distant(x$other, 0.9, 0, radius),
        `suspects-unrelated` = extract_distant(x$unrelated, 1, 0, radius))
    susp = susp[!sapply(susp, is.null)]
    if (length(susp) == 0)
        susp = NULL
    susp
}
