/* 
 * Vidjil browser, main configuration file
 * This file must be named 'js/conf.js' to be taken into account
 * */

var config = {

    /* Used for the 'align' script 
     * If this is not defined, the 'align' button is not available
     */    
    "cgi_address" : "default", // Public test server
    // "cgi_address" : "http://127.0.1.1/cgi-bin/",

    /* The following options control how the user may have access to .vidjil files.
     * Any combination of 1), 2) and 3) should work
     */

    /* 1) Patient database */
    "use_database" : true,
    "db_address" : "default",

    /* 2) and 3) Static files
    /* - relative paths if Vidjil browser is online on the same server
     * - absolute paths to a CORS active server only if browser is offline or on another server 
     */

    /* 2) Menu with a list of static files */
    "file_menu" : {
        "path" : "/data/",
        "file" : [
            "Stanford_S22.fasta"
        ]
    },
    
    /* 3) Static file autoload, possibly with an .analysis file */
    "autoload" : "data/Stanford-S22.vidjil"
    // "autoload_analysis" : "data/Stanford-S22.analysis"
}
