
σpr             = 0.25
Rcaps = [0.72, 1.1, 1.2, 1.3, 1.4]
# Rcaps = [2.28, 2.33, 2.47, 2.59, 2.67]

prior1           = OrderedDict{String, Distribution}()
prior1["γ"]      = Normal(mostdln(72.0, σpr*72.0)[1], mostdln(72.0, σpr*72.0)[2])
prior1["Rap"]    = Normal(mostdln(Rcaps[1], σpr*Rcaps[1])[1], mostdln(Rcaps[1], σpr*Rcaps[1])[2])

prior2           = OrderedDict{String, Distribution}()
prior2["γ"]      = Normal(mostdln(72.0, σpr*72.0)[1], mostdln(72.0, σpr*72.0)[2])
prior2["Rap"]    = Normal(mostdln(Rcaps[2], σpr*Rcaps[2])[1], mostdln(Rcaps[2], σpr*Rcaps[2])[2])

prior3           = OrderedDict{String, Distribution}()
prior3["γ"]      = Normal(mostdln(72.0, σpr*72.0)[1], mostdln(72.0, σpr*72.0)[2])
prior3["Rap"]    = Normal(mostdln(Rcaps[3], σpr*Rcaps[3])[1], mostdln(Rcaps[3], σpr*Rcaps[3])[2])

prior4           = OrderedDict{String, Distribution}()
prior4["γ"]      = Normal(mostdln(72.0, σpr*72.0)[1], mostdln(72.0, σpr*72.0)[2])
prior4["Rap"]    = Normal(mostdln(Rcaps[4], σpr*Rcaps[4])[1], mostdln(Rcaps[4], σpr*Rcaps[4])[2])

prior5           = OrderedDict{String, Distribution}()
prior5["γ"]      = Normal(mostdln(72.0, σpr*72.0)[1], mostdln(72.0, σpr*72.0)[2])
prior5["Rap"]    = Normal(mostdln(Rcaps[5], σpr*Rcaps[5])[1], mostdln(Rcaps[5], σpr*Rcaps[5])[2])