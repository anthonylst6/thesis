java -Dmf.cfg=config.cfg -jar aterm.jar nogui "asset.query :action get-values :where namespace>='%~1' :xpath -ename namespace namespace :xpath -ename name name :xpath -ename csum content/csum :xpath -ename size content/size :xpath -ename content content/type :xpath -ename ctime ctime :size infinity :output-format csv :out file:checksums.csv"
