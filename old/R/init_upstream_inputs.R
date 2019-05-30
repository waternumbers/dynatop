# configure any flows external to the channel
init_upstream_inputs <- function(upstream_inputs, groups, dt)
{
  upstream_inputs <- lapply(upstream_inputs,
    function(upstream_input)
    {
      # number of time steps before the flow reaches the outlet (dt in hours)
      idelay <- round(upstream_input$flow_dist_to_outlet/max(groups$vchan, na.rm=TRUE)/dt)

      nobs <- nrow(upstream_input$flow)
      iend <- nobs - idelay
      # don't worry too much about what happens in the first few time steps; replicate the first input
      # (or could make zero)
      upstream_input$qshift <- rbind(upstream_input$flow[rep(1,idelay),],
                      upstream_input$flow[1:iend,])

      return(upstream_input)

    }
  )
  return(upstream_inputs)
}
