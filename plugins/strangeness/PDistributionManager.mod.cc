//this line is added to the PDistributionManager.cc

PStrangenessPlugin *strange = new PStrangenessPlugin("strangeness","Strangeness extension");
AddPlugin(strange);
Enable("strangeness");
PluginInfo("Baryonic resonances with strangeness available");

