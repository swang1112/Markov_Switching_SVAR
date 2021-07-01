require(magrittr)
data      = read.table('data.txt') %>% ts(start = c(1973, 2), end = c(2007, 12), frequency = 12)
var3      = vars::VAR(data, p = 24)
data_time = time(data)[-(1:24)]
Bcoef     = vars::Bcoef(var3)[,1:(var3$K*var3$p)]
u         = var3 %>% resid