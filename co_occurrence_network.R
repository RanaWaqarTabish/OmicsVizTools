# Load necessary libraries
library(igraph)          # For network analysis
library(ggraph)          # For network visualization
library(tidygraph)       # For network data manipulation
library(dplyr)           # For data manipulation
library(psych)           # For correlation calculations

# ------------------------------
# Step 1: Load your dataset
# ------------------------------
# Example: Replace with your actual data loading step
# The data should be in a matrix or data frame format where rows are samples and columns are features (e.g., taxa, genes).
data <- read.csv("your_data.csv", row.names = 1)  # Replace 'your_data.csv' with your file name

# ------------------------------
# Step 2: Calculate Spearman correlations
# ------------------------------
# Compute Spearman correlations between all features
cor_matrix <- corr.test(data, method = "spearman", adjust = "BH")

# Extract correlation values and p-values
correlations <- cor_matrix$r
p_values <- cor_matrix$p

# ------------------------------
# Step 3: Apply thresholds for significance
# ------------------------------
# Thresholds: p-value < 0.05 and absolute correlation > 0.2
correlations[p_values > 0.05 | abs(correlations) < 0.2] <- 0

# ------------------------------
# Step 4: Create and simplify the network
# ------------------------------
# Build a network from the correlation matrix
network <- graph_from_adjacency_matrix(
  correlations, mode = "undirected", weighted = TRUE, diag = FALSE
)

# Remove isolated nodes (nodes with no connections)
network <- delete_vertices(network, V(network)[degree(network) == 0])

# ------------------------------
# Step 5: Analyze and extract key network properties
# ------------------------------
# Convert the network to a tidygraph object
network_tidy <- as_tbl_graph(network) %>%
  mutate(
    degree = centrality_degree(),
    community = as.factor(group_louvain())  # Community detection using Louvain algorithm
  )

# Identify the top nodes (hubs) in each community based on degree
top_nodes_per_community <- network_tidy %>%
  activate(nodes) %>%
  as_tibble() %>%
  group_by(community) %>%
  arrange(desc(degree)) %>%
  slice_head(n = 3) %>%  # Select top 3 nodes per community
  pull(name)

# ------------------------------
# Step 6: Visualize the network
# ------------------------------
# Set a random seed for reproducibility
set.seed(1234)

# Customize colors for communities (adjust as needed)
num_communities <- length(unique(network_tidy %>% activate(nodes) %>% pull(community)))
community_colors <- RColorBrewer::brewer.pal(min(12, num_communities), "Set3")

# Create the network plot
p <- ggraph(network_tidy, layout = "fr") +
  # Edges: Width and color based on correlation strength
  geom_edge_link(
    aes(width = abs(weight), color = weight),
    alpha = 0.6,
    show.legend = TRUE
  ) +
  scale_edge_color_gradient2(
    low = "#D73027", mid = "white", high = "#1A9850", midpoint = 0,
    name = "Correlation"
  ) +
  
  # Nodes: Size based on degree, color based on community
  geom_node_point(
    aes(size = degree, fill = community),
    shape = 21, color = "black", stroke = 0.3, alpha = 0.9
  ) +
  
  # Node labels: Display only top nodes from each community
  geom_node_text(
    aes(label = ifelse(name %in% top_nodes_per_community, name, "")),
    repel = TRUE, fontface = "bold", size = 3, color = "navy"
  ) +
  
  # Customize plot elements
  labs(
    title = "Feature Co-occurrence Network",
    subtitle = "Louvain Clustering | Spearman Correlations",
    caption = "Edge color: Correlation strength | Node size: Degree | Node color: Community"
  ) +
  
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.grid = element_blank(),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# ------------------------------
# Step 7: Save the plot
# ------------------------------
ggsave("Network.png", p, width = 12, height = 8, dpi = 300)
ggsave("Network.svg", p, width = 12, height = 8)
