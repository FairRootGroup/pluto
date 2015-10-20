//TITLE Access and print basic particle properties

{

//Make sure that data base is filled:
makeStaticData();

//The use of the wrapper functions is always
//a good starting point:
listParticle(-1);
//-1 means print all

listParticle(17);
//print out eta values with pid=17

//as all these are wrapper functions,
//one can also access the data base: 
makeDataBase()->ListEntries(-1,1,"*name,pid,*width");
//again: -1 means all, option=1 means one-liner
//the string also acts as a filter:
//only entries where at least one of
//the params are existing, are printed
//the names with "*" are NOT taken into
//account when filtering but printed
//in addition

//Print out which parameters are available:
makeDataBase()->Print();

}
