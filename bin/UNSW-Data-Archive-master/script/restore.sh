#!/usr/bin/env bash

java -Dmf.cfg=config.cfg -cp aterm.jar arc.mf.command.Execute "asset.preparation.request.create :migrate online :where namespace>='$1' :send-notification-email true :asset-service -name asset.content.validate"
