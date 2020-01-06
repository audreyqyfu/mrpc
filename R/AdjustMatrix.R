
# Adjust the columns of the input matrix to have

# the same ordering as in the reference matrix

AdjustMatrix <- function (reference, input) {
  
  if (any (colnames (reference) != colnames (input))) {
    
    Order_node <- match (colnames (reference), colnames (input))
    
    return (input[Order_node, Order_node])
    
  } else {
    
    return (input)
  }
}