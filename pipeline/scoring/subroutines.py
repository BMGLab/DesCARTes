

# Function to calculate CPS for one weight combination
def calculate_cps(row, weights):
    alpha, beta, gamma, delta, lambda_, eta, theta, mu, nu = weights
    AvgLogFC_Eo = row["AvgLogFC_Eo"]
    AvgLogFC_Ee = row["AvgLogFC_Ee"]
    B = row["B"]
    C = row["n_upregulated_norm"]
    L = row["L"]
    CV_intra = row["CV_intra"]
    CV_inter = row["CV_inter"]
    U_intra = row["U_intra"]
    U_inter = row["U_inter"]

    score = (alpha * AvgLogFC_Eo +
             beta * AvgLogFC_Ee +
             gamma * C +
             delta * L -
             lambda_ * np.log2(B + 1) -
             eta * CV_intra -
             theta * CV_inter +
             mu * U_intra +
             nu * U_inter)
    return score


def calculate_gene_score(df, weights):
    alpha, beta, gamma, delta, lambda_, eta, theta, mu, nu = weights
    B_safe = df["B"].clip(lower=1e-6)
    df["Score"] = (
        alpha * df["AvgLogFC_Eo"] +
        beta * df["AvgLogFC_Ee"] +
        gamma * df["n_upregulated_norm"] +
        delta * df["L"] -
        lambda_ * np.log2(B_safe + 1) -
        eta * df["CV_intra"] -
        theta * df["CV_inter"] +
        mu * df["U_intra"] +
        nu * df["U_inter"]
    )
    return df
