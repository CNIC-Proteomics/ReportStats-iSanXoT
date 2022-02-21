library("shiny")

args = commandArgs(trailingOnly=TRUE)
SRC_HOME = args[1]
setwd(SRC_HOME)

# Browser config
chrome.portable = file.path(getwd(),
                            "ChromiumPortable/chrome.exe")
launch.browser = function(appUrl, browser.path=chrome.portable) {
  message('Browser path: ', browser.path)
  shell(sprintf('"%s" --app=%s', browser.path, appUrl), wait=FALSE)
}

runApp("ShinyApp", launch.browser = launch.browser)