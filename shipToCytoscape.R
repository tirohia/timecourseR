install.packages("RJSONIO")
install.packages("httr")

library(RJSONIO)
library(igraph)
library(httr)

# Basic settings
port.number = 1234
base.url = paste("http://localhost:", toString(port.number), "/v1", sep="")

# Load utility functions
source("../code/R/toCytoscape.R")

# 0. Make sure Cytoscape REST API module is running
version.url = paste(base.url, "version", sep="/")
cytoscape.version = GET(version.url)
cy.version = fromJSON(rawToChar(cytoscape.version$content))
print(cy.version)

# 1. Create simple directed graph with Barabasi-Albert model
graph1 <- barabasi.game(200)
par(mfrow=c(1,1))
plot(graph1)

# 2. Calculate some statistics and assign then to the graph
graph1$name = "Scale-Free Network (BA Model)"
graph1$density = graph.density(graph1)

V(graph1)$degree <- degree(graph1)
V(graph1)$closeness <- closeness(graph1)
V(graph1)$betweenness <- betweenness(graph1)
V(graph1)$page_rank <- page.rank(graph1)$vector
V(graph1)$community <- label.propagation.community(graph1)$membership

E(graph1)$betweenness <- edge.betweenness(graph1)

graph1

# 3. Build a custom Visual Style programatically
style.name = "MyCustomStyle3"

# Defaults
def.node.color <- list(
  visualProperty = "NODE_FILL_COLOR",
  value = "#00aabb"
)

def.node.border.width <- list(
  visualProperty = "NODE_BORDER_WIDTH",
  value = 0
)

def.edge.target.arrow <- list(
  visualProperty="EDGE_TARGET_ARROW_SHAPE",
  value="ARROW"
)

defaults <- list(def.node.color, def.node.border.width,def.edge.target.arrow)

# Visual Mappings
min.betweenness = min(V(graph1)$betweenness)
max.betweenness = max(V(graph1)$betweenness)

mappings = list()

point1 = list(
  value=min.betweenness,
  lesser= "20.0",
  equal="20.0",
  greater="20.0"
)

point2 = list(
  value=max.betweenness,
  lesser="200.0",
  equal="200.0",
  greater="200.0"
)

node.size.continuous.points = list(point1, point2)

node.size = list(
  mappingType="continuous",
  mappingColumn="betweenness",
  mappingColumnType="Integer",
  visualProperty="NODE_SIZE",
  points = node.size.continuous.points
)

node.label = list(
  mappingType="passthrough",
  mappingColumn="name",
  mappingColumnType="String",
  visualProperty="NODE_LABEL"
)

mappings = list(node.size, node.label)

style <- list(title=style.name, defaults = defaults, mappings = mappings)
style.JSON <- toJSON(style)

style.url = paste(base.url, "styles", sep="/")
POST(url=style.url, body=style.JSON, encode = "json")

# Convert to Cytoscape style JSON object
cygraph <- toCytoscape(graph)

# Send it to Cytoscape!
network.url = paste(base.url, "networks", sep="/")
res <- POST(url=network.url, body=cygraph, encode="json")
network.suid = unname(fromJSON(rawToChar(res$content)))


# 5. Apply layouts and Visual Style
apply.layout.url = paste(base.url, "apply/layouts/force-directed", toString(network.suid), sep="/")
apply.style.url = paste(base.url, "apply/styles", style.name, toString(network.suid), sep="/")

res <- GET(apply.layout.url)
res <- GET(apply.style.url)
