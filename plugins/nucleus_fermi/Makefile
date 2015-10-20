###

LOCAL_CLASSES    := PFermiDistributions PFermiMomentumGA PNucleusFermiPlugin 
LOCAL_COLLECTION := PDistributionManager.mod.cc



### do not change something below
#LOCAL_PATH :=$(strip $(subst Makefile, , $(lastword $(MAKEFILE_LIST))))
LOCAL_PATH :=$(strip $(subst Makefile, , $(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST))))
PLUGIN_CLASSES := $(PLUGIN_CLASSES) $(addprefix $(LOCAL_PATH),$(LOCAL_CLASSES))
PLUGIN_CLASSES_NAMES := $(PLUGIN_CLASSES_NAMES) $(LOCAL_CLASSES)
PLUGIN_COLLECTION := $(PLUGIN_COLLECTION) $(addprefix $(LOCAL_PATH),$(LOCAL_COLLECTION))
