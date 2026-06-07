# Archive ABC Filtering

This workflow performs archive-based sequential ABC filtering over saved broad
`+/-5%` ODE archive rows.

It lowers distance tolerances across rounds, but it does not generate new ODE
proposals. It is therefore an archive-based ABC approximation rather than a
full adaptive ABC-SMC sampler. Only ODE-confirmed archive rows are interpreted
scientifically.

